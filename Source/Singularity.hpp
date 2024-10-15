#pragma once

#include "Assertion.hpp"

template <class T>
class Singularity {
protected:
  static T* singularity_;

public:
  Singularity() {
    ASSERTION(!singularity_, "Singularity instance has already been created!");
    singularity_ = static_cast<T*>(this);
  }

  virtual ~Singularity() {
    ASSERTION(singularity_, "Singularity instance has to be created first!");
    singularity_ = nullptr;
  }

  /** Note: Scott Meyers mentions in his Effective Modern
   *        C++ book, that deleted functions should generally
   *        be public as it results in better error messages
   *        due to the compilers behavior to check accessibility
   *        before deleted status.
   * https://stackoverflow.com/questions/1008019/c-singleton-design-pattern
   */
  Singularity(const Singularity& other)                  = delete;
  const Singularity& operator=(const Singularity& other) = delete;

  static T& getSingularity() {
    ASSERTION(singularity_, "Singularity instance has to be created first!");
    return *singularity_;
  }

  static T* getSingularityPtr() { return singularity_; }
};

template <class T>
T* Singularity<T>::singularity_ = nullptr;
