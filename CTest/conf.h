#pragma once
#include "Testdata_path.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

// static const std::string model_name("/beam");
static const int num_model_cells = 865;

template<class F, typename... Types>
double
time_check(F&& f, Types&&... args)
{
  auto start = std::chrono::system_clock::now();
  f(std::forward<Types>(args)...);
  auto end = std::chrono::system_clock::now();
  double dura = std::chrono::duration<double, std::milli>(end - start).count();
  return dura;
}

template<typename T>
double
random_real(T lb, T rb)
{
  std::random_device rd{};
  std::default_random_engine re(rd());
  return std::uniform_real_distribution<T>(lb, rb)(re);
}