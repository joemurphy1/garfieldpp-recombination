#ifndef G_EXCEPTION
#define G_EXCEPTION

#include <exception>
#include <string>

namespace Garfield {

class Exception : public std::exception {
 public:
  explicit Exception(const std::string& message) : m_message(message) {}
  explicit Exception(const char* const message) : m_message(message) {}
  const char* what() const noexcept override { return m_message.c_str(); }
  Exception(Exception&& other) = delete;
  Exception& operator=(Exception&& other) = delete;
  Exception& operator=(const Exception& other) = delete;
  Exception(const Exception& other) = delete;

 private:
  std::string m_message;
};

}  // namespace Garfield

#endif