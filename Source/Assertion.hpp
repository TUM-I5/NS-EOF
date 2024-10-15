#pragma once

struct Assertion {
  static bool throwAssertionException;

  typedef bool (*Handler)();
  static Handler handler;

  class AssertionException: public std::exception {
    const std::string message_;

  public:
    AssertionException(const std::string& message);
    virtual const char* what() const noexcept override;
  };

  struct ScopedAssertionExceptionEnabler {
    const bool originalValue;

    ScopedAssertionExceptionEnabler();
    ~ScopedAssertionExceptionEnabler();
  };

  template <class T0, class... TArgs>
  static void stringify(std::stringstream& stream, T0&& t0, TArgs&&... rest) {
    stream << t0;
    if constexpr (sizeof...(TArgs) > 0) {
      stream << ", ";
    }
    stringify(stream, std::forward<TArgs>(rest)...);
  }

  static void stringify(std::stringstream&) {}

  template <class... TArgs>
  static void fire(const char* message, const char* file, const char* func, const int line, TArgs&&... args) {
    std::stringstream stream{};
    stringify(stream, std::forward<TArgs>(args)...);
    fireParams(message, file, func, line, stream.str().data());
  }

  static void fireParams(
    const char* message, const char* file, const char* func, const int line, const char* params = ""
  );
};

// Helper cast, performing a static_cast, but checking that the cast is valid using dynamic_cast in assert and debug
// builds.
template <class TDerived, class TBase>
static inline TDerived assertionCast(TBase* value) {
  static_assert(std::is_pointer<TDerived>::value, "Must be a pointer type");
  ASSERTION(!value || dynamic_cast<TDerived>(value) != nullptr, value);
  return static_cast<TDerived>(value);
}

#ifndef NDEBUG
// https://docs.microsoft.com/en-us/cpp/intrinsics/debugbreak?view=vs-2017
#define ASSERTION(x, ...) \
  if (!(x)) { \
    Assertion::fire(#x, __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    assert(true); \
  }

#define ASSERTION_EQUALS(lhs, rhs, ...) \
  if ((lhs) != (rhs)) { \
    std::stringstream ss; \
    ss << #lhs << "==" << #rhs; \
    Assertion::fire(ss.str().data(), __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    assert(false); \
  }
#else
#define ASSERTION(x, ...)
#define ASSERTION_EQUALS(lhs, rhs, ...)
#endif
