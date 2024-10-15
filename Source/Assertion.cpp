#include "StdAfx.hpp"

#include "Assertion.hpp"

bool Assertion::throwAssertionException = true;

Assertion::Handler Assertion::handler = nullptr;

Assertion::AssertionException::AssertionException(const std::string& message):
  message_(message) {}

const char* Assertion::AssertionException::what() const noexcept { return message_.data(); }

Assertion::ScopedAssertionExceptionEnabler::ScopedAssertionExceptionEnabler():
  originalValue(throwAssertionException) {

  throwAssertionException = true;
}

Assertion::ScopedAssertionExceptionEnabler::~ScopedAssertionExceptionEnabler() {
  throwAssertionException = originalValue;
}

static std::vector<std::string> getStackTrace() {
#ifndef _MSC_VER
  constexpr std::uint32_t TRACE_SIZE = 10;
  void*                   buffer[TRACE_SIZE]{};
  const int               nptrs = backtrace(buffer, TRACE_SIZE);

  char** strings = backtrace_symbols(buffer, nptrs);
  if (!strings) {
    return {"No backtrace could be generated!"};
  }

  std::vector<std::string> trace;
  for (int i = 0; i < nptrs; ++i) {
    trace.push_back(strings[i]);
  }

  free(strings);
  return trace;

#else
  return {"No backtrace could be generated!"};
#endif
}

void Assertion::fireParams(
  const char* message, const char* file, const char* func, const int line, const char* params
) {
  static std::mutex            mutex;
  std::unique_lock<std::mutex> lock(mutex);

  std::stringstream ss;
  ss << "[" << file << ":" << line << "]: Assertion \"" << message << "\" failed in function \"" << func << "()\""
     << "\n";
  if (strlen(params) != 0) {
    ss << "Assertion parameters: " << params << "\n";
  }
  ss << "Stack trace: ";
  const auto& trace = getStackTrace();
  for (const auto& s : trace) {
    ss << s;
  }

  if (throwAssertionException) {
    throw AssertionException(ss.str());
  } else {
#ifdef _MSC_VER
    if (IsDebuggerPresent()) {
      __debugbreak(); // Make the MSVC++ Debugger stop here if it is attached
    }
#else
    raise(SIGTRAP);
#endif
  }
}
