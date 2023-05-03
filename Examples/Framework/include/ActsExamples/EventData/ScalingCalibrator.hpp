// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// TODO move this somewhere else?

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <ActsExamples/EventData/MeasurementCalibration.hpp>

#include <TFile.h>
#include <TH2D.h>

namespace ActsExamples {

class ScalingCalibrator : public MeasurementCalibrator {
 public:
  struct ConstantTuple {
    double x_offset;
    double x_scale;
    double y_offset;
    double y_scale;
  };

  struct MapTuple {
    TH2D x_offset;
    TH2D x_scale;
    TH2D y_offset;
    TH2D y_scale;

    ConstantTuple at(size_t sizeLoc0, size_t sizeLoc1) const {
      ConstantTuple ct;
      ct.x_offset =
          x_offset.GetBinContent(x_offset.FindFixBin(sizeLoc0, sizeLoc1));
      ct.x_scale =
          x_scale.GetBinContent(x_scale.FindFixBin(sizeLoc0, sizeLoc1));
      ct.y_offset =
          y_offset.GetBinContent(y_offset.FindFixBin(sizeLoc0, sizeLoc1));
      ct.y_scale =
          y_scale.GetBinContent(y_scale.FindFixBin(sizeLoc0, sizeLoc1));
      return ct;
    }
  };

  ScalingCalibrator(const char* path);

  void calibrate(
      const MeasurementContainer& measurements,
      const ClusterContainer* clusters, const Acts::GeometryContext& /*gctx*/,
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
          trackState) const;

  bool needsClusters() { return true; }

 private:
  std::map<Acts::GeometryIdentifier, MapTuple> m_calib_maps;
  std::bitset<3> m_mask;
};

}  // namespace ActsExamples
