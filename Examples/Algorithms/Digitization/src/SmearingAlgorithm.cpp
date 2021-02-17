// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <algorithm>
#include <stdexcept>
#include <type_traits>

ActsExamples::SmearingAlgorithm::SmearingAlgorithm(
    ActsExamples::SmearingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SmearingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements output collection");
  }
  if (m_cfg.outputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links output collection");
  }
  if (m_cfg.outputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-particles map output collection");
  }
  if (m_cfg.outputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map output collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (not m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }
  if (m_cfg.smearers.empty()) {
    throw std::invalid_argument("Missing smearers configuration");
  }
  // create the smearers from the configuration
  std::vector<std::pair<Acts::GeometryIdentifier, Smearer>> smearersInput;
  for (size_t i = 0; i < m_cfg.smearers.size(); ++i) {
    Acts::GeometryIdentifier geoId = m_cfg.smearers.idAt(i);
    // copy so we can sort in-place
    SmearingConfig geoCfg = m_cfg.smearers.valueAt(i);

    // make sure the configured input parameter indices are sorted and unique
    std::sort(geoCfg.begin(), geoCfg.end(),
              [](const ParameterSmearingConfig& lhs,
                 const ParameterSmearingConfig& rhs) {
                return lhs.index < rhs.index;
              });
    auto dup = std::adjacent_find(geoCfg.begin(), geoCfg.end(),
                                  [](const ParameterSmearingConfig& lhs,
                                     const ParameterSmearingConfig& rhs) {
                                    return lhs.index == rhs.index;
                                  });
    if (dup != geoCfg.end()) {
      std::invalid_argument(
          "Smearer configuration contains duplicate parameter indices");
    }

    switch (geoCfg.size()) {
      case 1u:
        smearersInput.emplace_back(geoId, makeSmearer<1u>(geoCfg));
        break;
      case 2u:
        smearersInput.emplace_back(geoId, makeSmearer<2u>(geoCfg));
        break;
      case 3u:
        smearersInput.emplace_back(geoId, makeSmearer<3u>(geoCfg));
        break;
      case 4u:
        smearersInput.emplace_back(geoId, makeSmearer<4u>(geoCfg));
        break;
      default:
        throw std::invalid_argument("Unsupported smearer size");
    }
  }
  m_smearers = Acts::GeometryHierarchyMap<Smearer>(std::move(smearersInput));
}

ActsExamples::ProcessCode ActsExamples::SmearingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // retrieve input
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);

  // prepare output containers
  IndexSourceLinkContainer sourceLinks;
  MeasurementContainer measurements;
  IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
  IndexMultimap<Index> hitSimHitsMap;
  sourceLinks.reserve(simHits.size());
  measurements.reserve(simHits.size());
  hitParticlesMap.reserve(simHits.size());
  hitSimHitsMap.reserve(simHits.size());

  // setup random number generator
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  for (auto simHitsGroup : groupByModule(simHits)) {
    // manual pair unpacking instead of using
    //   auto [moduleGeoId, moduleSimHits] : ...
    // otherwise clang on macos complains that it is unable to capture the local
    // binding in the lambda used for visiting the smearer below.
    Acts::GeometryIdentifier moduleGeoId = simHitsGroup.first;
    const auto& moduleSimHits = simHitsGroup.second;

    auto smearerItr = m_smearers.find(moduleGeoId);
    if (smearerItr == m_smearers.end()) {
      ACTS_DEBUG("No smearing function present for module " << moduleGeoId);
      continue;
    }
    const Acts::Surface* surfacePtr =
        m_cfg.trackingGeometry->findSurface(moduleGeoId);
    if (not surfacePtr) {
      // this is either an invalid geometry id or a misconfigured smearer
      // setup; both cases can not be handled and should be fatal.
      ACTS_ERROR("Could not find surface " << moduleGeoId
                                           << " for configured smearer");
      return ProcessCode::ABORT;
    }

 
    // run the smearer. iterate over the hits for this surface inside the
    // visitor so we do not need to lookup the variant object per-hit.
    std::visit(
        [&](const auto& smearer) {
          using ThisSmearer = std::decay_t<decltype(smearer)>;
          using ThisMeasurement =
              Acts::Measurement<IndexSourceLink, Acts::BoundIndices,
                                ThisSmearer::size()>;

	  // Intermediary storage
	  std::vector<ThisMeasurement> unmerged_measurements;
	  std::vector<std::vector<ActsFatras::Hit>> associated_simhits; // FIXME

          for (auto h = moduleSimHits.begin(); h != moduleSimHits.end(); ++h) {
            const auto& simHit = *h;
            // const auto simHitIdx = simHits.index_of(h);

            auto res = smearer(rng, simHit, *surfacePtr, ctx.geoContext);
            if (not res.ok()) {
              // ignore un-smearable measurements
              // TODO log this or at least count invalid hits?
              return;
            }
            const auto& [par, cov] = res.value();

	    Index hitIdx = unmerged_measurements.size();
	    IndexSourceLink sourceLink(moduleGeoId, hitIdx);
            ThisMeasurement meas(sourceLink, smearer.indices, par, cov);

	    unmerged_measurements.emplace_back(std::move(meas));
	    std::vector<ActsFatras::Hit> v{simHit};
	    associated_simhits.push_back(v);
	  }

	  // FIXME: this should be optimized
	  for (size_t i = 0; i < unmerged_measurements.size(); i++) {
	    if (associated_simhits.at(i).empty()) // already merged
	      continue;

	    // N.B.
	    // This method only works if the covariance is diagonal

	    auto& fst = unmerged_measurements.at(i);
	    auto& par_fst = fst.parameters();
	    auto& covm_fst = fst.covariance();
	    Acts::Vector2 pos_fst { par_fst[Acts::ePos0], par_fst[Acts::ePos1] };
	    Acts::SymMatrix2 cov_fst = Acts::SymMatrix2::Zero();
	    cov_fst(0, 0) = covm_fst(Acts::ePos0, Acts::ePos0);
	    cov_fst(1, 1) = covm_fst(Acts::ePos1, Acts::ePos1);

	    for (size_t j = 0; j < unmerged_measurements.size(); j++) {
	      if (i == j or associated_simhits.at(j).empty()) // already merged
		continue;

	      auto& snd = unmerged_measurements.at(j);
	      auto& par_snd = snd.parameters();
	      auto& covm_snd = snd.covariance();
	      Acts::Vector2 pos_snd { par_snd[Acts::ePos0], par_snd[Acts::ePos1] };
	      Acts::SymMatrix2 cov_snd = Acts::SymMatrix2::Zero();
	      cov_snd(0, 0) = covm_snd(Acts::ePos0, Acts::ePos0);
	      cov_snd(1, 1) = covm_snd(Acts::ePos1, Acts::ePos1);

	      auto diff = (pos_snd - pos_fst).normalized();

	      auto proj_pos_fst = pos_fst.dot(diff);
	      auto proj_pos_snd = pos_fst.dot(diff);

	      double proj_cov_fst = diff.transpose() * cov_fst * diff;
	      double proj_cov_snd = diff.transpose() * cov_snd * diff;

	      // TODO: configurable criterion
	      if ((proj_pos_fst + std::sqrt(proj_cov_fst)) > (proj_pos_snd - std::sqrt(proj_cov_snd))) {
		// MERGE
		for (auto &sh : associated_simhits.at(j))
		  associated_simhits.at(i).push_back(std::move(sh));
		associated_simhits.at(j).clear();

		// update position
		// FIXME doens't change the measurement
		pos_fst += pos_snd;
		pos_fst *= 0.5;
		cov_fst += cov_snd;
		cov_fst *= 0.5;
	      }
	    }
	  }

	  for (size_t i = 0; i < unmerged_measurements.size(); i++) {
	    if (associated_simhits.at(i).empty())
	      continue; // Merged into another
	    Index hitIdx = measurements.size();
	    IndexSourceLink sourceLink(moduleGeoId, hitIdx);
	    sourceLinks.emplace_hint(sourceLinks.end(), std::move(sourceLink));
	    // FIXME: the source link will be wrong here
	    measurements.emplace_back(std::move(unmerged_measurements[i]));
	    for (auto &sh : associated_simhits.at(i)) {
	      // FIXME: does this work like I think it does?
	      // const auto simHitIdx = simHits.index_of(sh); // FIXME BROKEN
	      hitParticlesMap.emplace_hint(hitParticlesMap.end(), hitIdx, sh.particleId());
	      // hitSimHitsMap.emplace_hint(hitSimHitsMap.end(), hitIdx, simHitIdx); // FIXME BROKEN
	    }
	  }
        },
        *smearerItr);
  }

  ctx.eventStore.add(m_cfg.outputSourceLinks, std::move(sourceLinks));
  ctx.eventStore.add(m_cfg.outputMeasurements, std::move(measurements));
  ctx.eventStore.add(m_cfg.outputMeasurementParticlesMap,
                     std::move(hitParticlesMap));
  ctx.eventStore.add(m_cfg.outputMeasurementSimHitsMap,
                     std::move(hitSimHitsMap));
  return ProcessCode::SUCCESS;
}
