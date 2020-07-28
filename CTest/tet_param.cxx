#include "LatticeCore/Analysis/field/frame.h"
#include "LatticeCore/Mesh/geometry.h"
#include "conf.h"
#include "gtest/gtest.h"
using namespace std;

TEST(tet_param_test, test_one_tet)
{
  V3d a({ 0, 0, 0 }), b({ 1, 0, 0 }), c({ 0, -1, 0 }), d({ 0, 0, -2 });
  LatticeCore::Tetrahedron tet = { a, b, c, d };
  Matrix_3 frame =
    LatticeCore::QuaternionFrame(Eigen::Quaterniond::UnitRandom())
      .to_rotation_matrix();
  V3d rho{ 1.0, 1.0, 1.0 };
  V4d f(0, 1, -1, -2);
  auto der = tet.get_grad_operator();
  V3d x1 = der * f;
  cout << x1 << endl;

  auto inv_jac = tet.get_inv_jac_mat();
  V3d x2 = inv_jac * V3d(f(1) - f(0), f(2) - f(0), f(3) - f(0));
  EXPECT_NEAR((x1 - x2).norm(), 0, 1e-4);
}