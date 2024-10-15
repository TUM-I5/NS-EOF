#pragma once

template <class T>
class Singleton {
  static bool       isInstantiated_;
  static std::mutex instanceMutex_;

protected:
  explicit Singleton() = default;
  virtual ~Singleton() = default;

public:
  Singleton(const Immovable&)            = delete;
  Singleton& operator=(const Singleton&) = delete;

  static const std::shared_ptr<T>& getInstance() {
    // This implementation is thread safe
    // https://stackoverflow.com/questions/1661529/is-meyers-implementation-of-the-singleton-pattern-thread-safe
    static_assert(not std::is_same<T, nullptr_t>::value);
    static std::shared_ptr<T> object(new T());
    isInstantiated_ = true;
    return object;
  }

  static bool isInstantiated() { return isInstantiated_; }
};

template <class T>
bool Singleton<T>::isInstantiated_ = false;

template <class T>
std::mutex Singleton<T>::instanceMutex_;

#define BEFRIEND_SINGLETON \
  template <class> \
  friend class Singleton
