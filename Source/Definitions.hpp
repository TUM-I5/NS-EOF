#pragma once

#ifndef ENABLE_PETSC
#define PETSC_COMM_WORLD MPI_COMM_WORLD
#endif

// Datatype for the type of data stored in the structures
#ifdef ENABLE_SINGLE_PRECISION
using RealType = float;
#define MY_MPI_FLOAT MPI_FLOAT
#else
using RealType = double;
#define MY_MPI_FLOAT MPI_DOUBLE
#endif

static constexpr RealType MY_FLOAT_MAX = std::numeric_limits<RealType>::max();
static constexpr RealType MY_FLOAT_MIN = std::numeric_limits<RealType>::min();

#define LIKELY(x) __builtin_expect(bool(x), 1)
#define UNLIKELY(x) __builtin_expect(bool(x), 0)
#define NORETURN [[noreturn]]
#define BIT(x) (1 << x)
#define DEPRECATED __attribute__((deprecated))
#define NOINLINE __attribute__((noinline))
#define FALLTHROUGH [[fallthrough]]
#define _STR(x) #x
#define TO_STRING(x) _STR(x)

#ifdef _MSC_VER
#ifndef FORCEINLINE
#define FORCEINLINE inline __forceinline
#endif
#define TALIGN(n) __declspec(align(n)) // alignas(n)
#else
#define FORCEINLINE inline __attribute__((always_inline))
#define ALIGN(n) __attribute__((aligned(n)))
#endif

// A fallback implementation of uintptr_t for systems that lack it
struct fallback_uintptr {
  unsigned char value[sizeof(void*)];
};
