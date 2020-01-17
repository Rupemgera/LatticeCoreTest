#include "LatticeCore/Implementation/field/frame_field.h"
#include "LatticeCore/Implementation/mesh/geometry.h"
#include "conf.h"
#include "gtest/gtest.h"
#include <cmath>

TEST(geometry_test, tet_test) {
  V3d a({0, 0, 0}), b({0, 0, 1}), c({1, 0, 0}), d({0, 1, 0});
  LatticeCore::Tetrahedron tet = {a, b, c, d};

  // std::cerr << "coeffs : \n" << tet.get_coeffs() << std::endl;

  auto test_interpolation = [&tet](V3d &&coor) {
    // assume f(x,y,z) = x+y-z, testing whether interpolation is correct.
    V4d phi(0, -1, 1, 1);
    auto coeffs = tet.interpolate(coor);
    EXPECT_DOUBLE_EQ(coor(0) + coor(1) - coor(2), coeffs.dot(phi)) << coeffs;
    // test grad
    auto gopt = tet.get_grad_operator();
    V3d grad = gopt * phi;
    EXPECT_DOUBLE_EQ(0.0, (grad - V3d(1, 1, -1)).squaredNorm()) << grad;
  };

  // test a point inside tet
  test_interpolation(V3d(0.25, 0.25, 0.25));

  // test vertices
  test_interpolation(V3d(0, 1, 0));
}

TEST(parameterization_test, test_one_tet_para) {
  // auto parmtor = LatticeCore::Parameterizor();
  V3d a({0, 0, 0}), b({1, 0, 0}), c({0, -1, 0}), d({0, 0, -2});
  LatticeCore::Tetrahedron tet = {a, b, c, d};
  Matrix_3 frame;
  frame << 1, 0, 0, 0, 1, 0, 0, 0, 1;

  Eigen::Matrix<double, 4, 3> phi;
  for (size_t i = 0; i < 3; i++) {
    phi.col(i) = tet.parameterize(frame, i);
  }
  auto tframe = tet.grad_vector(phi);
  EXPECT_NEAR(LatticeCore::Frame3d::metric(tframe, frame), 0.0, 1e-7)
      << phi << std::endl
      << tframe << std::endl
      << frame << std::endl;
}

TEST(parameterization_test, test_one_tet_para2) {
  // auto parmtor = LatticeCore::Parameterizor();
  V3d a({0, 8, 0}), b({8, 0, 8}), c({0, 0, 0}), d({0, 3.12246, 0.25932});
  LatticeCore::Tetrahedron tet = {a, b, c, d};
  Matrix_3 frame;
  frame << 0.0284146, -0.992387, 0.119839, 0.26573, 0.123073, 0.956159,
      -0.963629, 0.00467603, 0.2672040;

  Eigen::Matrix<double, 4, 3> phi;
  for (size_t i = 0; i < 3; i++) {
    phi.col(i) = tet.parameterize(frame, i);
  }

  Matrix_3 m = tet.grad_vector(phi);
  // std::cout << phi << std::endl;
  EXPECT_NEAR(LatticeCore::Frame3d::metric(m, frame), 0.0, 1e-7);
}

// given some values of nodes, parameterize tet
TEST(parameterization_test, test_fixed_para) {
  V3d a({0, 0, 0}), b({1, 1, 0}), c({-1, 1, 0}), d({0, 1, -2});
  LatticeCore::Tetrahedron tet = {a, b, c, d};
  Matrix_3 frame;
  frame << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  std::vector<size_t> ids = {0, 1};
  std::vector<double> vals = {1.0, 2.0};
  auto phi = tet.parameterize(frame, ids, vals);
  // std::cout << "fixed_para: phi: \n" << phi << std::endl;
  EXPECT_DOUBLE_EQ(phi(1) - phi(0), 1.0) << phi << std::endl;
  EXPECT_DOUBLE_EQ(phi(2) - phi(0), -1.0) << phi << std::endl;
  EXPECT_DOUBLE_EQ(phi(3) - phi(0), 0.0) << phi << std::endl;
}

TEST(parameterization_test, test_fixed_para_axis_y) {
  V3d a({0, 0, 0}), b({1, 0, 0}), c({0, -1, 0}), d({0, 0, -2});
  LatticeCore::Tetrahedron tet = {a, b, c, d};
  Matrix_3 frame;
  frame << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  std::vector<size_t> ids = {0, 2};
  std::vector<double> vals[] = {{1.0, 1.0}, {1.0, 0.0}, {0.0, 0.0}};
  Eigen::Matrix<double, 4, 3> phi;
  for (size_t i = 0; i < 3; i++) {
    phi.col(i) = tet.parameterize(frame, ids, vals[i], i);
  }
  auto tframe = tet.grad_vector(phi);
  EXPECT_NEAR(LatticeCore::Frame3d::metric(tframe, frame), 0.0, 1e-7)
      << tframe << std::endl
      << frame << std::endl;
}

// test whether Frame3d::metric is transform-invariant
TEST(frame_test, test_frame_metric) {
  auto seed_frame =
      LatticeCore::QuaternionFrame(Eigen::Quaterniond::UnitRandom());
  auto frame1 = seed_frame.to_rotation_matrix();
  Matrix_3 frame2 = frame1 * -1;
  frame2.col(2) = frame1.col(2);
  EXPECT_NEAR(LatticeCore::Frame3d::metric(frame1, frame2), 0.0, 1e-9);
}

TEST(frame_test, test_metric_range) {
  Matrix_3 frame1 = Matrix_3::Identity();
  Matrix_3 frame2 = frame1;
  const int n_case = 20;
  double step = PI / n_case;
  // std::vector<double> diffs(n_case);
  double maxd = -1e20;
  for (size_t i = 0; i < n_case; i++) {
    double angle = step * i;
    Eigen::AngleAxisd Rz(angle, V3d::UnitZ());
    Eigen::AngleAxisd Rx(angle, V3d::UnitX());
    frame2 = Rz * Rx * frame2;
    double d = LatticeCore::Frame3d::metric(frame1, frame2);
    // diffs[i] = d;
    // std::cout << angle << " : " << d << std::endl;
    if (d > maxd)
      maxd = d;
  }
  EXPECT_GT(maxd, 1) << maxd;
}
