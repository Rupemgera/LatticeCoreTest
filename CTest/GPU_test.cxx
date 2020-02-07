#include "LatticeCore/Parellel/cuda_export.h"
#include "conf.h"
#include "gtest/gtest.h"

TEST(gpu_test, device_info_test) {
  int i = device_info();
  EXPECT_EQ(i, 0);
}

TEST(gpu_test, vector_additon_test) {
  constexpr int size = 1000;
  float a[size];
  float b[size];
  float gpu_res, cpu_res;
  for (int i = 0; i < size; ++i) {
    a[i] = 1.0;
    b[i] = 1.0;
  }
  auto gpu_add = [a, b, &gpu_res]() { gpu_res = vector_add(a, b, size); };
  auto cpu_add = [a, b, &cpu_res]() {
    cpu_res = 0.0;
    for (int i = 0; i < size; ++i) {
      cpu_res += a[i] * b[i];
    }
  };
  std::cout << "GPU time: " << time_check(gpu_add) << " ms\n";
  std::cout << "CPU time: " << time_check(cpu_add) << " ms\n";
  EXPECT_NEAR(gpu_res, cpu_res, 1e-9);
}