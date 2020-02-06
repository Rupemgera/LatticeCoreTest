#include "LatticeCore/Implementation/field/frame_field.h"
#include "LatticeCore/Implementation/mesh/geometry.h"
#include "conf.h"
#include "dlib/optimization.h"
#include "gtest/gtest.h"
#include <cmath>

TEST(parameterization_test, test_optimization_one_tet) {
  V3d a({0, 0, 0}), b({1, 0, 0}), c({0, -1, 0}), d({0, 0, -2});
  LatticeCore::Tetrahedron tet = {a, b, c, d};
  using Vector = dlib::matrix<double, 0, 1>;
  Matrix_3 frame =
      LatticeCore::QuaternionFrame(Eigen::Quaterniond::UnitRandom())
          .to_rotation_matrix();
  V3d rho{1.0, 1.0, 1.0};
  auto f_obj = [&tet, &frame, &rho](const Vector &m) {
    V4d x{m(0), m(1), m(2), m(3)};
    return tet.objective(frame, x, rho, 0);
  };

  auto f_grad = [&tet, &frame, &rho](const Vector &m) {
    V4d x{m(0), m(1), m(2), m(3)};
    auto dx = tet.objective_gradient(frame, x, rho, 0);
    Vector r(4);
    for (size_t i = 0; i < 4; i++) {
      r(i) = dx(i);
    }
    return r;
  };

  Vector x(4);
  for (size_t i = 0; i < 4; i++) {
    x(i) = i;
  }

  dlib::find_min(dlib::lbfgs_search_strategy(40),
                 dlib::objective_delta_stop_strategy(1e-7).be_verbose(), f_obj,
                 f_grad, x, -1);

  V4d y;
  for (size_t i = 0; i < 4; i++) {
    y(i) = x;
  }

  V3d v = tet.grad_operator() * y;

  EXPECT_NEAR((v.dot(frame.col(0))) * (v.dot(frame.col(0))), 1.0, 1e-9)
      << v << std::endl
      << "y: " << y << std::endl
      << frame.col(0) << std::endl;
}

TEST(optimization_test, test_local_minima) {
  using Vector = dlib::matrix<double, 0, 1>;
  auto f = [](const Vector &m) {
    double s = 0.0;
    for (size_t i = 0; i < 2; i++) {
      s += m(i) * m(i);
    }
    s -= 1.0;
    return s * s + m(0) * m(0) * m(1) * m(1);
  };
  auto df = [](const Vector &m) {
    double s = 0.0;
    for (size_t i = 0; i < 2; i++) {
      s += m(i) * m(i);
    }
    s -= 1.0;
    Vector res(2);
    res(0) = s * m(0) * 4.0 + m(0) * m(1) * m(1) * 2.0;
    res(1) = s * m(1) * 4.0 + m(0) * m(0) * m(1) * 2.0;

    return res;
  };
  Vector x(2);
  x(0) = 1;
  x(1) = 1;
  dlib::find_min(dlib::lbfgs_search_strategy(40),
                 dlib::objective_delta_stop_strategy(1e-7).be_verbose(), f, df,
                 x, -1);
  std::cout << x << std::endl;
  EXPECT_NEAR(x(0) * x(0) + x(1) * x(1), 1.0, 1e-9);
}

// TEST(optimization_test, test_gpu) {
// }