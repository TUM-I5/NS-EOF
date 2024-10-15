#pragma once

#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#endif

#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <assert.h>
#include <bitset>
#include <chrono>
#include <cmath>
#include <complex>
#include <csignal>
#include <cstring>
#include <deque>
#include <exception>
#include <filesystem>
#include <float.h>
#include <forward_list>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <mpi.h>
#include <mutex>
#include <numeric>
#include <random>
#include <set>
#include <setjmp.h>
#include <signal.h>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <thread>
#include <time.h>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>

// https://en.cppreference.com/w/cpp/feature_test
#if __has_include(<format>)
#include <format>
#endif

#ifdef _MSC_VER
#include <direct.h>
#include <process.h>
#include <Windows.h>
#include <Winsock2.h>
#else
#include <execinfo.h>
#include <unistd.h>
#endif

#ifdef ENABLE_PETSC
#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#endif
