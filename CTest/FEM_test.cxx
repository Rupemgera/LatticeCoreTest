#include "LatticeCore/Analysis/FEM/truss.h"
#include "conf.h"
#include <gtest/gtest.h>

using namespace LatticeCore;

TEST(fem_test, test_bar)
{
  std::vector<V3d> p = {
    V3d(0, 0, 0), V3d(400, 0, 0), V3d(400, 300, 0), V3d(0, 300, 0)
  };
  std::vector<std::array<int, 2>> connect = {
    { 0, 1 }, { 1, 2 }, { 0, 2 }, { 2, 3 }
  };
  std::vector<double> A = { 100, 100, 100, 100 };
  std::vector<std::array<double, 3>> f = {
    { 0, 0, 0 }, { 20000, 0, 0 }, { 0, -25000, 0 }, { 0, 0, 0 }
  };
  std::vector<std::array<int, 2>> fixed = { { 0, 0 }, { 0, 1 }, { 0, 2 },
                                            { 1, 1 }, { 1, 2 }, { 2, 2 },
                                            { 3, 0 }, { 3, 1 }, { 3, 2 } };
  double E = 29.5e4;

  TrussFEM solver(p, connect, A, E);
  auto q = solver.analyse_displacement_vector(f, fixed);
  EXPECT_EQ(q.size(), 12);
  EXPECT_NEAR(q[3], 0.2712, 1e-4);
  EXPECT_NEAR(q[6], 0.0565, 1e-4);
  EXPECT_NEAR(q[7], -0.2225, 1e-4);
}