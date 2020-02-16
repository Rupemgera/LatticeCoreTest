#include "LatticeCore/Parellel/cuda_export.h"
#include "LatticeCore/Parellel/cuda_optimization.h"
#include "conf.h"
#include "gtest/gtest.h"

TEST(gpu_test, device_info_test) {
  int i = device_info(true);
  EXPECT_EQ(i, 0);
}

TEST(gpu_test, rosenbrock_test) {
  constexpr int N = 4006000;
  // constexpr int N = (1 << 22) + (1 << 10);
  float *x = new float[N];
  for (long long i = 0; i < N; ++i) {
    x[i] = 0;
  }

  double c_sum = 0.0;
  auto start = std::chrono::system_clock::now();
  for (long long id = 0; id < N - 1; ++id) {
    c_sum += (x[id + 1] - x[id] * x[id]) * (x[id + 1] - x[id] * x[id]) * 100.0 +
             (x[id] - 1.0) * (x[id] - 1.0);
    // if (id % 1000000 == 0)
    //   std::cout << "id: " << id << " sum: " << std::setprecision(12) << c_sum
    //             << std::endl;
  }
  auto end = std::chrono::system_clock::now();
  double dura = std::chrono::duration<double, std::milli>(end - start).count();
  std::cout << "CPU time: " << dura << "ms\n";

  double g_sum = 0.0;
  auto rosen_gpu = [x, &g_sum]() { g_sum = Rosenbrock_GPU(x, N); };

  std::cout << "GPU time: " << time_check(rosen_gpu) << "ms\n";
  delete[] x;

  EXPECT_NEAR(g_sum, c_sum, 1e-9);
}