// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DIGITIZATION_PLANARMODULECLUSTER_H
#define ACTS_DIGITIZATION_PLANARMODULECLUSTER_H 1

#include "ACTS/Digitization/DigitizationCell.hpp"
#include "ACTS/EventData/Measurement.hpp"
#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Identifier.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {

template <ParID_t... params>
using Measurement_t = Measurement<Identifier, params...>;

class PlanarModuleCluster : public Measurement_t<ParDef::eLOC_0, ParDef::eLOC_1>
{
public:
  /// Constructor from DigitizationCells
  ///
  /// @param mSurface is the module surface
  /// @param cIdentifier is the channel identifier of the local position
  /// @param cov is the covariance matrix
  /// @param loc0 is the local position in the first coordinate
  /// @param loc1 is the local position in the second coordinate
  /// @param dCells is the vector of digitization cells
  /// - optional truth information
  /// @param barcodes particle barcodes of simulated particles
  /// @param tparameters particle impact position, momentum
  PlanarModuleCluster(const Surface&                mSurface,
                      const Identifier&             cIdentifier,
                      ActsSymMatrixD<2>             cov,
                      double                        loc0,
                      double                        loc1,
                      std::vector<DigitizationCell> dCells,
                      std::vector<barcode_type>     barcodes    = {},
                      std::vector<Vector3D>         tparameters = {})
    : Measurement_t<ParDef::eLOC_0, ParDef::eLOC_1>(mSurface,
                                                    cIdentifier,
                                                    std::move(cov),
                                                    loc0,
                                                    loc1)
    , m_digitizationCells(dCells)
    , m_barcodes(barcodes)
    , m_tparameters(tparameters)
  {
  }

  /// access to the digitization cells
  ///
  /// @return the vector to the digitization cells
  const std::vector<DigitizationCell>&
  digitizationCells() const;

  /// access to the contributing barcodes
  ///
  /// @return the vector of the particle barcode
  const std::vector<barcode_type>&
  barcodes() const;

  /// access to the truth parameters (if present)
  ///
  /// @return a vector of truth parameters at the cluster
  const std::vector<Vector3D>&
  truthParameters() const;

private:
  std::vector<DigitizationCell> m_digitizationCells;  /// the digitization cells
  std::vector<barcode_type>     m_barcodes;           /// barcodes of particles
  std::vector<Vector3D>         m_tparameters;        /// truth parameters
};

inline const std::vector<DigitizationCell>&
PlanarModuleCluster::digitizationCells() const
{
  return m_digitizationCells;
}

inline const std::vector<barcode_type>&
PlanarModuleCluster::barcodes() const
{
  return m_barcodes;
}

inline const std::vector<Vector3D>&
PlanarModuleCluster::truthParameters() const
{
  return m_tparameters;
}
}

#endif  // ACTS_DIGITIZATION_PLANARMODULECLUSTER_H
