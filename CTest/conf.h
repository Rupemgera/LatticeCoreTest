#pragma once
#include "Testdata_path.h"
#include <chrono>
#include <fstream>
#include <iostream>

static const std::string model_name("/beam");
static const int num_model_cells = 865;

template <class F> double time_check(F &f) {
  auto start = std::chrono::system_clock::now();
  f();
  auto end = std::chrono::system_clock::now();
  double dura = std::chrono::duration<double, std::milli>(end - start).count();
  return dura;
}