#include "LatticeCore/Component/optimization/PSO.hpp"
#include "conf.h"
#include "gtest/gtest.h"

using namespace PSO;
TEST(PSO_test, Rosenbrock)
{
  auto obj = [](const Eigen::VectorXd& x) {
    auto n = x.rows();
    double sum = 0.0;
    for (int i = 0; i < n - 1; ++i) {
      sum += (x(i + 1) - x(i) * x(i)) * (x(i + 1) - x(i) * x(i)) * 100.0 +
             (x(i) - 1.0) * (x(i) - 1.0);
    }
    return sum;
  };

  auto constraint = [](const Eigen::VectorXd& x) { return true; };

  constexpr int P = 40;
  constexpr int N = 20;
  ParticalSwarm particals(P, Eigen::VectorXd(N));
  ParticalSwarm velocity(P, Eigen::VectorXd(N));
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < N; j++) {
      particals[i](j) = random_real<double>(0, 1);
      velocity[i](j) = random_real<double>(0, 1);
    }
  }

  int max_iter = 200;
  LinearRatio w(0.8, 0.4, max_iter);
  Ratio c1(2.0), c2(2.0);
  ParticalSwarmOptimization pso(&w, &c1, &c2);
  auto x = pso.optimize(particals, velocity, obj, constraint, N, max_iter);
  EXPECT_NEAR(x.sum(), N, 1e-3 * N);
  EXPECT_NEAR(obj(x), 0.0, 1e-3);
  for (int i = 0; i < N; ++i) {
    if (std::abs(x(i) - 1.0) > 0.5) {
      std::cout << i << " " << x(i) << std::endl;
    }
  }
}

TEST(PSO_test, Quad)
{
  auto obj = [](const Eigen::VectorXd& x) {
    auto n = x.rows();
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      sum += x(i) * x(i);
    }
    return sum;
  };

  auto constraint = [](const Eigen::VectorXd& x) { return true; };

  constexpr int P = 40;
  constexpr int N = 20;
  ParticalSwarm particals(P, Eigen::VectorXd(N));
  ParticalSwarm velocity(P, Eigen::VectorXd(N));
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < N; j++) {
      particals[i](j) = random_real<double>(-1, 1);
      velocity[i](j) = random_real<double>(-1, 1);
    }
  }

  int max_iter = 200;
  LinearRatio w(0.8, 0.4, max_iter);
  Ratio c1(2.0), c2(2.0);
  ParticalSwarmOptimization pso(&w, &c1, &c2);
  auto x = pso.optimize(particals, velocity, obj, constraint, N, max_iter);
  EXPECT_NEAR(x.sum(), 0.0, 1e-3 * N);
  EXPECT_NEAR(obj(x), 0.0, 1e-3);
}