
#include "LatticeCore/Analysis/field/frame_field.h"
#include "LatticeCore/Mesh/geometry.h"
#include "conf.h"
#include "modules/quaternion/quat.h"
#include "gtest/gtest.h"
#include <cmath>
#include <tuple>

using QFrame = LatticeCore::QuaternionFrame;
using namespace std;

// test whether Frame3d::metric is transform-invariant
TEST(frame_test, test_frame_metric)
{
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

TEST(frame_test, test_metric_range)
{
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

class QuaternionParameterTestClass
  : public ::testing::TestWithParam<std::tuple<int, int>>
{
  // : public ::testing::TestWithParam<int> {

protected:
  int i1, i2;
  virtual void SetUp()
  {
    i1 = std::get<0>(GetParam());
    i2 = std::get<1>(GetParam());
  }
  virtual void TearDown() {}
};

TEST_P(QuaternionParameterTestClass, Gao_quaternion_tes)
{
  double qs[7][4] = { { -1.330099286227237332e-01,
                        -9.265940911951956460e-01,
                        3.358196234978461092e-02,
                        3.501485411288546290e-01 },
                      { 1.829444081596818053e-02,
                        5.325000829396204782e-01,
                        8.144714121275935248e-02,
                        -8.423035903359292753e-01 },
                      { -3.516813571819379924e-02,
                        4.202375266151780475e-01,
                        1.612714680540463519e-01,
                        -8.922752585643946022e-01 },
                      { -6.995875727057111471e-01,
                        1.242956833085924057e-01,
                        6.891627718730093388e-01,
                        1.420650734378147151e-01 },
                      { -4.629436714740193803e-01,
                        5.067084649481913283e-01,
                        5.880723770538690554e-01,
                        -4.279025215375747981e-01 },
                      { -5.724395513527590768e-02,
                        -7.892805381874992143e-01,
                        -6.788004529779008249e-03,
                        -6.113209342345683472e-01 },
                      { -1.395174782939428426e-01,
                        -4.909076455163013941e-01,
                        1.643267347745258899e-01,
                        8.441216032435628902e-01 } };

  // int i1 = 5;
  // int i2 = GetParam();

  QFrame qf1(qs[i1][1], qs[i1][2], qs[i1][3], qs[i1][0]);
  QFrame qf2(qs[i2][1], qs[i2][2], qs[i2][3], qs[i2][0]);

  auto rqf = QFrame::best_match(qf2, qf1);
  auto v1 = rqf.coeffs();
  for (int i = 0; i < 4; ++i)
    cout << v1(i) << ' ';
  cout << endl;

  Quaternion q1(qs[i1][1], qs[i1][2], qs[i1][3], qs[i1][0]);
  Quaternion q2(qs[i2][1], qs[i2][2], qs[i2][3], qs[i2][0]);
  auto rq = Quaternion::findRotation(q1, q2);
  auto rq2 = Quaternion::applyRotation(q1, q2);
  Eigen::Vector4d v2 = rq;
  Eigen::Vector4d v3 = rq2;
  Eigen::Vector4d v4 = (qf1 * rqf).coeffs();
  for (int i = 0; i < 4; ++i)
    cout << v2(i) << ' ';
  cout << endl;

  double d1 = (v1 + v2).head<3>().norm();
  double d2 = (v1 - v2).head<3>().norm();
  double d3 = v1(3) - v2(3) + std::min(d1, d2);
  EXPECT_NEAR(d3, 0, 1e-7);
  EXPECT_NEAR((v3 - v4).norm(), 0, 1e-7);
}

// INSTANTIATE_TEST_CASE_P(inn, QuaternionParameterTestClass,
//                         ::testing::Range(0, 6));
INSTANTIATE_TEST_SUITE_P(pairs,
                         QuaternionParameterTestClass,
                         ::testing::Combine(::testing::Range(0, 7),
                                            ::testing::Range(0, 7)));