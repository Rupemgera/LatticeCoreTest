#include "LatticeCore/lattice_core.h"
#include "Testdata_path.h"
#include <string>

using namespace LatticeCore;
using namespace std;

class LatticeConsole {
public:
  LatticeWrapper lattice;

  LatticeConsole() {
    string model_name = "beam";
    string path = string(TEST_DATA_PATH) + "/";
    string mesh_file = path + model_name + ".ovm";
    string stress_file = path + model_name + ".csv";
    string frame_file = path + model_name + ".txt";
    lattice.read_mesh(mesh_file);
    lattice.read_stress_field(stress_file);
  }

  void smooth() { lattice.request_smoothed_stress_field(1.0, 4.0); }
};

int main() {
  LatticeConsole console;
  // console.smooth();
  return 0;
}