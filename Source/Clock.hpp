#pragma once

class Clock {
private:
  typedef std::chrono::high_resolution_clock              ClockType;
  typedef std::chrono::duration<std::uint64_t, std::nano> DurationType;

  ClockType::time_point start_;

public:
  Clock();
  ~Clock() = default;

  std::uint64_t getTime() const;
  std::string   getDate() const;

  static void        sleep(std::uint32_t msec);
  static std::string formatHMS(std::int64_t t);
};
