#include "LatticeCore/Parellel/cuda_export.h"
#include "LatticeCore/Parellel/cuda_optimization.h"
#include "conf.h"
#include "gtest/gtest.h"

TEST(gpu_test, device_info_test) {
  int i = device_info();
  EXPECT_EQ(i, 0);
}

TEST(gpu_test, rosenbrock_test) {
  constexpr int N = 8000000;
  float *x = new float[N];
  for (int i = 0; i < N; ++i) {
    x[i] = 0.0;
  }

  float c_sum = 0.0;
  auto rosen_cpu = [x, &c_sum]() {
    for (int id = 0; id < N - 1; ++id) {
      c_sum +=
          (x[id + 1] - x[id] * x[id]) * (x[id + 1] - x[id] * x[id]) * 100.0 +
          (x[id] - 1.0) * (x[id] - 1.0);
    }
  };

  std::cout << "CPU time: " << time_check(rosen_cpu) << "ms\n";

  float g_sum = 0.0;
  auto rosen_gpu = [x, &g_sum]() { g_sum = Rosenbrock_GPU(x, N); };

  std::cout << "GPU time: " << time_check(rosen_gpu) << "ms\n";
  delete[] x;

  EXPECT_NEAR(g_sum, c_sum, 1e-9);
}