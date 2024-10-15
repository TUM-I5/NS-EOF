#include "StdAfx.hpp"

#include "Clock.hpp"

Clock::Clock():
  start_(Clock::ClockType::now()) {}

std::uint64_t Clock::getTime() const {
  Clock::DurationType elapsed = Clock::ClockType::now() - start_;
  return elapsed.count();
}

std::string Clock::getDate() const {
  auto        now = std::chrono::system_clock::now();
  std::time_t t   = std::chrono::system_clock::to_time_t(now);
  std::string str = std::ctime(&t);
  str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
  return str;
}

void Clock::sleep(std::uint32_t msec) {
#ifdef _MSC_VER
  std::this_thread::sleep_for(std::chrono::microseconds(msec * 1000));
#else
  usleep(msec * 1000);
#endif
}

std::string Clock::formatHMS(std::int64_t t) {
  int hours = int(t / (std::int64_t(3600) * std::int64_t(1000000000)));
  t -= hours * std::int64_t(3600) * int64_t(1000000000);
  int minutes = int(t / (std::int64_t(60) * std::int64_t(1000000000)));
  t -= minutes * std::int64_t(60) * int64_t(1000000000);
  int seconds = int(t / std::int64_t(1000000000));

  std::ostringstream str;
  str << std::setfill('0') << std::setw(2) << hours << ":" << std::setw(2) << minutes << ":" << std::setw(2) << seconds;
  return str.str();
}
