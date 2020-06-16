#include "LatticeCore/Component/optimization/PSO.hpp"
#include "conf.h"
#include "gtest/gtest.h"

using namespace PSO;

TEST(PSO_test, dice)
{
  auto rd_engine = std::default_random_engine{ std::random_device{}() };
  std::uniform_real_distribution<double> urd(0.0, 1.0);
  auto dice = std::bind(urd, rd_engine);
  constexpr int N = 1000;
  std::vector<double> digs;
  digs.reserve(N);
  double sum = 0.0;
  for (int i = 0; i < N; ++i) {
    digs.push_back(dice());
    sum += digs.back();
  }
  double variance = 0.0;
  double average = sum / N;
  for (int i = 0; i < N; ++i) {
    variance += (digs[i] - average) * (digs[i] - average) / N;
  }
  std::sort(digs.begin(), digs.end());
  EXPECT_GE(digs.front(), 0.0);
  EXPECT_LE(digs.back(), 1.0);
  EXPECT_LE(digs.front() - digs.back(), 1.0);
  EXPECT_NEAR(average, 0.5, 0.1);
  EXPECT_NEAR(variance, 1.0 / 12, 0.01);
  std::cout << "average: " << average << std::endl;
  std::cout << "variance: " << variance << std::endl;
}

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

  constexpr int P = 50;
  constexpr int N = 20;
  ParticalSwarm particals(P, Eigen::VectorXd(N));
  ParticalSwarm velocity(P, Eigen::VectorXd(N));
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < N; j++) {
      particals[i](j) = random_real<double>(0, 2);
      velocity[i](j) = random_real<double>(0, 2);
    }
  }

  int max_iter = 200;
  LinearRatio w(0.9, 0.6, max_iter);
  Ratio c1(2.0), c2(2.0);
  ParticalSwarmOptimization pso(&w, &c1, &c2);
  auto x = pso.optimize(particals, velocity, obj, constraint, N, max_iter);
  EXPECT_NEAR(x.sum(), N, 1e-3 * N);
  EXPECT_NEAR(obj(x), 0.0, 1e-3);
  // for (int i = 0; i < N; ++i) {
  //   if (std::abs(x(i) - 1.0) > 0.5) {
  //     std::cout << i << " " << x(i) << std::endl;
  //   }
  // }
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

  constexpr int P = 50;
  constexpr int N = 60;
  ParticalSwarm particals(P, Eigen::VectorXd(N));
  ParticalSwarm velocity(P, Eigen::VectorXd(N));
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < N; j++) {
      particals[i](j) = random_real<double>(-1, 1);
      velocity[i](j) = random_real<double>(-1, 1);
    }
  }

  int max_iter = 200;
  LinearRatio w(0.9, 0.6, max_iter);
  Ratio c1(2.5), c2(2.5);
  ParticalSwarmOptimization pso(&w, &c1, &c2);
  auto x = pso.optimize(particals, velocity, obj, constraint, N, max_iter);
  EXPECT_NEAR(x.sum(), 0.0, 1e-3 * N);
  EXPECT_NEAR(obj(x), 0.0, 1e-3);
}