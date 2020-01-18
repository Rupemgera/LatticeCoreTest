
#include "LatticeCore/Implementation/field/frame_field.h"
#include "LatticeCore/Implementation/mesh/geometry.h"
#include "conf.h"
#include "gtest/gtest.h"
#include <cmath>

// test whether Frame3d::metric is transform-invariant
TEST(frame_test, test_frame_metric) {
  auto seed_frame =
      LatticeCore::QuaternionFrame(Eigen::Quaterniond::UnitRandom());
  auto frame1 = seed_frame.to_rotation_matrix();
  Matrix_3 frame2 = frame1 * -1.0;
  frame2.col(2) = frame1.col(2);
  Matrix_3 frame3;
  frame3.col(0) = frame1.col(1);
  frame3.col(1) = frame1.col(2);
  frame3.col(2) = frame1.col(0) * -1.0;
  EXPECT_NEAR(LatticeCore::Frame3d::metric(frame1, frame2), 0.0, 1e-9);
  EXPECT_NEAR(LatticeCore::Frame3d::metric(frame1, frame3), 0.0, 1e-9)
      << frame1 << std::endl
      << frame3 << std::endl;
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

TEST(frame_test, test_objective) {}
