#include "LatticeCore/DataInterface/mesh_interface.hpp"
#include "LatticeCore/lattice_core.h"
#include "conf.h"
#include <gtest/gtest.h>

class MeshTest : public ::testing::Test {
public:
  const std::string model_name = "beam";

protected:
  MeshTest() : lattice() {
    path = std::string(TEST_DATA_PATH) + "/";
    std::string mesh_file = path + model_name + ".ovm";
    std::string stress_file = path + model_name + ".csv";
    std::string frame_file = path + model_name + ".txt";
    lattice.read_mesh(mesh_file);
    lattice.read_stress_field(stress_file);
    // lattice.read_frames(frame_file);
    mesh_data = lattice.request_mesh_data();
  }

  virtual ~MeshTest() {}

  virtual void SetUp() {}

  virtual void TearDown() {}

  LatticeCore::LatticeWrapper lattice;
  LatticeCore::MeshInterface mesh_data;
  std::string path;
};

TEST_F(MeshTest, Readin) {
  auto centers = mesh_data.request_cell_centers();
  size_t n = centers.size();

  EXPECT_EQ(n, num_model_cells);
}

TEST_F(MeshTest, frame_optimzaiton) {
  lattice.request_smoothed_stress_field(1.0, 4.0);

  EXPECT_EQ(1, 1);
}
// TEST_F(MeshTest, GPU) {
//   int i = lattice.GPU_info();
//   EXPECT_EQ(i, 0);
// }