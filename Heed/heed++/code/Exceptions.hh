#ifndef HEED_EXCEPTIONS
#define HEED_EXCEPTIONS

#include <stdexcept>
#include <string>
#if __cplusplus >= 202002L || _MSVC_LANG >= 202002L
#include <source_location>
#endif

namespace Heed {

class Exception : public std::runtime_error {
 public:
#if __cplusplus >= 202002L || _MSVC_LANG >= 202002L
  explicit Exception(
      const char* const message,
      const std::source_location& location = std::source_location::current())
      : std::runtime_error(message), m_source_location(location) {}
  virtual const char* what() const noexcept {
    try {
      static thread_local std::string full_message;
      full_message.clear();
      full_message = std::runtime_error::what();
      full_message += " in " + std::string(m_source_location.function_name()) +
                      " at " + std::string(m_source_location.file_name()) +
                      ":" + std::to_string(m_source_location.line()) + ":" +
                      std::to_string(m_source_location.column());
      return full_message.c_str();
    } catch (...) {
      return std::runtime_error::what();
    }
  }
#else
  explicit Exception(const char* const message) : std::runtime_error(message) {}
#endif
 private:
#if __cplusplus >= 202002L || _MSVC_LANG >= 202002L
  std::source_location m_source_location;
#endif
};

}  // namespace Heed

#endif