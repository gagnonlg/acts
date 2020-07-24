// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This file contains systematic tests to verify the local->global->local
// transformation roundtrip for all available surfaces with a large range of
// possible parameters.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <cmath>
#include <limits>

namespace {

namespace bdata = boost::unit_test::data;
using namespace Acts;

constexpr auto eps = 8 * std::numeric_limits<double>::epsilon();
const GeometryContext geoCtx;

void runTest(const Surface& surface, double l0, double l1, double phi,
             double theta) {
  const Vector3D dir = makeDirectionUnitFromPhiTheta(phi, theta);

  // convert local-to-global
  Vector3D sentinel = Vector3D::Random();
  Vector3D pos = sentinel;
  surface.localToGlobal(geoCtx, Vector2D(l0, l1), dir, pos);
  BOOST_TEST(pos != sentinel, "Position was not changed");
  BOOST_TEST(std::isfinite(pos[0]),
             "Position " << pos.transpose() << " contains non-finite entries");
  BOOST_TEST(std::isfinite(pos[1]),
             "Position " << pos.transpose() << " contains non-finite entries");
  BOOST_TEST(std::isfinite(pos[2]),
             "Position " << pos.transpose() << " contains non-finite entries");
  BOOST_TEST(surface.isOnSurface(geoCtx, pos, dir),
             "Position " << pos.transpose() << " is not on the surface");

  // convert global-to-local
  Vector2D loc = Vector2D::Zero();
  bool validTransform = surface.globalToLocal(geoCtx, pos, dir, loc);
  BOOST_TEST(validTransform);
  CHECK_CLOSE_OR_SMALL(loc[ePos0], l0, eps, eps);
  CHECK_CLOSE_OR_SMALL(loc[ePos1], l1, eps, eps);
}

std::shared_ptr<Transform3D> makeTransformIdentity() {
  return std::make_shared<Transform3D>(Transform3D::Identity());
}

// test datasets

// local positions
const auto posAngle = bdata::xrange(-M_PI, M_PI, 0.25);
const auto posPositiveNonzero = bdata::xrange(0.25, 1.0, 0.25);
const auto posPositive = bdata::make(0.0) + posPositiveNonzero;
const auto posSymmetric = bdata::xrange(-1.0, 1.0, 0.25);
// direction angles
const auto phis = bdata::xrange(-M_PI, M_PI, M_PI_4);
const auto thetasNoForwardBackward = bdata::xrange(M_PI_4, M_PI, M_PI_4);
const auto thetas = bdata::make({0.0, M_PI}) + thetasNoForwardBackward;

// different surfaces
// parameters must be chosen such that all possible local positions (as defined
// in the datasets above) represent valid points on the surface.
const auto cones = bdata::make({
    Surface::makeShared<ConeSurface>(makeTransformIdentity(),
                                     0.5 /* opening angle */),
});
const auto cylinders = bdata::make({
    Surface::makeShared<CylinderSurface>(makeTransformIdentity(),
                                         10.0 /* radius */, 100 /* half z */),
});
const auto discs = bdata::make({
    Surface::makeShared<DiscSurface>(makeTransformIdentity(),
                                     0 /* radius min */, 100 /* radius max */),
});
const auto perigees = bdata::make({
    Surface::makeShared<PerigeeSurface>(Vector3D(0, 0, -1.5)),
});
const auto planes = bdata::make({
    Surface::makeShared<PlaneSurface>(Vector3D(1, 2, 3), Vector3D::UnitX()),
    Surface::makeShared<PlaneSurface>(Vector3D(-2, -3, -4), Vector3D::UnitY()),
    Surface::makeShared<PlaneSurface>(Vector3D(3, -4, 5), Vector3D::UnitZ()),
});
const auto straws = bdata::make({
    Surface::makeShared<StrawSurface>(makeTransformIdentity(), 2.0 /* radius */,
                                      200.0 /* half z */),
});

}  // namespace

BOOST_AUTO_TEST_SUITE(SurfaceLocalToGlobalRoundtrip)

BOOST_DATA_TEST_CASE(ConeSurface,
                     cones* posAngle* posPositiveNonzero* phis* thetas, surface,
                     lphi, lz, phi, theta) {
  // TODO extend lz to zero after fixing the transform implementation
  // local parameter r*phi has limits that depend on the z position
  const auto r = lz * surface->bounds().tanAlpha();
  // local coordinates are singular at z = 0 -> normalize local phi
  runTest(*surface, (0 < lz) ? (r * lphi) : 0.0, lz, phi, theta);
}

BOOST_DATA_TEST_CASE(CylinderSurface,
                     cylinders* posSymmetric* posSymmetric* phis* thetas,
                     surface, lrphi, lz, phi, theta) {
  runTest(*surface, lrphi, lz, phi, theta);
}

BOOST_DATA_TEST_CASE(DiscSurface, discs* posPositive* posAngle* phis* thetas,
                     surface, lr, lphi, phi, theta) {
  // local coordinates are singular at r = 0 -> normalize local phi
  runTest(*surface, lr, (0 < lr) ? lphi : 0.0, phi, theta);
}

BOOST_DATA_TEST_CASE(
    PerigeeSurface,
    perigees* posSymmetric* posSymmetric* phis* thetasNoForwardBackward,
    surface, d0, z0, phi, theta) {
  // TODO extend theta to forward/back extreme cases fixing the transform
  runTest(*surface, d0, z0, phi, theta);
}

BOOST_DATA_TEST_CASE(PlaneSurface,
                     planes* posSymmetric* posSymmetric* phis* thetas, surface,
                     l0, l1, phi, theta) {
  runTest(*surface, l0, l1, phi, theta);
}

BOOST_DATA_TEST_CASE(
    StrawSurface,
    straws* posSymmetric* posSymmetric* phis* thetasNoForwardBackward, surface,
    lr, lz, phi, theta) {
  // TODO extend theta to forward/back extreme cases fixing the transform
  runTest(*surface, lr, lz, phi, theta);
}

BOOST_AUTO_TEST_SUITE_END()