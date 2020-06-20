#include "LatticeCore/Analysis/FEM/truss.h"
#include "conf.h"
#include <gmock/gmock.h>
// #include <gtest/gtest.h>

using namespace LatticeCore;
using namespace testing;

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

TEST(fem_test, test1)
{
  std::vector<V3d> p;
  for (int i = 0; i < 9; ++i) {
    p.push_back(V3d(i, (i % 2) * 2, 1));
  }
  std::vector<std::array<int, 2>> connect{ { 0, 1 }, { 0, 2 }, { 1, 2 },
                                           { 1, 3 }, { 2, 3 }, { 2, 4 },
                                           { 3, 4 }, { 3, 5 }, { 4, 5 },
                                           { 4, 6 }, { 5, 6 }, { 5, 7 },
                                           { 6, 7 }, { 6, 8 }, { 7, 8 } };
  double E = 3e7;
  std::vector<double> A;
  for (int i = 0; i < 15; ++i) {
    A.push_back(0.02 + 0.025 * (i % 2));
  }
  std::vector<std::array<int, 2>> fixed = { { 0, 0 }, { 0, 1 }, { 8, 1 } };
  for (int i = 0; i < 9; ++i) {
    fixed.push_back({ i, 2 });
  }
  std::vector<std::array<double, 3>> f(9, { 0, 0, 0 });
  f[1][0] = 15;
  f[2][1] = -5;
  f[3][1] = -7;
  f[6][1] = -10;
  TrussFEM solver(p, connect, A, E);
  auto q = solver.analyse_displacement_vector(f, fixed);
  double res[18];
  double decimal = 100000;
  for (int i = 0; i < 9; i++) {
    res[i * 2] = std::round(q[i * 3] * decimal) / decimal;
    res[i * 2 + 1] = std::round(q[i * 3 + 1] * decimal) / decimal;
  }
  double expect[] = { 0.00000, 0.00000,  0.00014, -0.00010, 0.00003, -0.00019,
                      0.00010, -0.00023, 0.00006, -0.00023, 0.00007, -0.00021,
                      0.00009, -0.00018, 0.00005, -0.00009, 0.00010, 0.00000 };
  EXPECT_THAT(res, Pointwise(DoubleEq(), expect));
}