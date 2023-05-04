// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include <ActsExamples/EventData/Measurement.hpp>

#include <cassert>

namespace ActsExamples {

/// Abstract base class for measurement-based calibration
class MeasurementCalibrator {
 public:
  virtual void calibrate(
      const MeasurementContainer& measurements,
      const ClusterContainer* clusters, const Acts::GeometryContext& gctx,
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
          trackState) const = 0;

  virtual ~MeasurementCalibrator() = default;
  virtual bool needsClusters() { return false; }
};

// Calibrator to convert an index source link to a measurement as-is
class PassThroughCalibrator : public MeasurementCalibrator {
 public:
  /// Find the measurement corresponding to the source link.
  ///
  /// @tparam parameters_t Track parameters type
  /// @param gctx The geometry context (unused)
  /// @param trackState The track state to calibrate
  void calibrate(
      const MeasurementContainer& measurements,
      const ClusterContainer* /*clusters*/,
      const Acts::GeometryContext& /*gctx*/,
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
          trackState) const;
};

// Executor class that pairs together a MeasurementCalibrator and a
// MeasurementContainer, to be passed to a Kalman Fitter, for instance
class MeasurementCalibratorAdapter {
 public:
  MeasurementCalibratorAdapter(const MeasurementCalibrator& calibrator,
                               const MeasurementContainer& measurements,
                               const ClusterContainer* clusters = nullptr);

  MeasurementCalibratorAdapter() = delete;

  void calibrate(
      const Acts::GeometryContext& gctx,
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy
          trackState) const;

 private:
  const MeasurementCalibrator& m_calibrator;
  const MeasurementContainer& m_measurements;
  const ClusterContainer* m_clusters;
};

}  // namespace ActsExamples
