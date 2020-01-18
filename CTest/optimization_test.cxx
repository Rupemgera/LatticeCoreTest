#include "LatticeCore/Implementation/field/frame_field.h"
#include "LatticeCore/Implementation/mesh/geometry.h"
#include "conf.h"
#include "gtest/gtest.h"
#include <cmath>

TEST(parameterization_test, test_optimization) {
  V3d a({0, 0, 0}), b({1, 0, 0}), c({0, -1, 0}), d({0, 0, -2});
  LatticeCore::Tetrahedron tet = {a, b, c, d};
}
