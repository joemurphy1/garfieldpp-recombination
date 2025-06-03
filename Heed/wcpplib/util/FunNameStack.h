#ifndef FUNNAMESTACK_H
#define FUNNAMESTACK_H

/*
Copyright (c) 1999 I. B. Smirnov

Permission to use, copy, modify, distribute and sell this file
and its documentation for any purpose is hereby granted without fee,
provided that the above copyright notice, this permission notice,
and notices about any modifications of the original text
appear in all copies and in supporting documentation.
It is provided "as is" without express or implied warranty.
*/

#include <cstdlib>
#include <iostream>
#include <list>

#include "wcpplib/stream/prstream.h"

// #define USE_BOOST_MULTITHREADING
// #define PRINT_MESSAGE_ABOUT_THREAD_INITIALIZATION
//  Prints a long message at initialization of FunNameStack.
//  The message is quite long and it is unlikely that
//  it can be visually missed even if it is interrupted by play of threads.
//  This message is not necessary in routine work.

/*
// Here there is a good place to switch off the
// initialization of the function names in all the programs
#ifdef FUNNAMESTACK
#undef FUNNAMESTACK
#endif
*/

// Switch on/off initialization of function names.
#define mfunname(string)

#define check_econd1(condition, a1, stream) \
  if (condition) {                          \
    stream << "ERROR:\n";                   \
    stream << '\n' << #condition << '\n';   \
    stream << #a1 << '=' << (a1) << '\n';   \
    spexit(stream);                         \
  }
#define check_wcond1(condition, a1, stream) \
  if (condition) {                          \
    stream << "WARNING:\n";                 \
    stream << '\n' << #condition << '\n';   \
    stream << #a1 << '=' << (a1) << '\n';   \
  }

#define check_econd11(a, signb, stream)     \
  if (a signb) {                            \
    stream << "ERROR:\n";                   \
    stream << '\n' << #a << #signb << '\n'; \
    stream << #a << '=' << (a) << '\n';     \
    spexit(stream);                         \
  }

#define check_econd12(a, sign, b, stream)                          \
  if (a sign b) {                                                  \
    stream << "ERROR:\n";                                          \
    stream << '\n' << #a << #sign << #b << '\n';                   \
    stream << #a << '=' << (a) << ' ' << #b << '=' << (b) << '\n'; \
    spexit(stream);                                                \
  }

// condition + additional any commands
#define check_econd11a(a, signb, add, stream) \
  if (a signb) {                              \
    stream << "ERROR:\n";                     \
    stream << '\n' << #a << #signb << '\n';   \
    stream << #a << '=' << (a) << '\n';       \
    stream << add;                            \
    spexit(stream);                           \
  }

#define check_econd12a(a, sign, b, add, stream)                    \
  if (a sign b) {                                                  \
    stream << "ERROR:\n";                                          \
    stream << '\n' << #a << #sign << #b << '\n';                   \
    stream << #a << '=' << (a) << ' ' << #b << '=' << (b) << '\n'; \
    stream << add;                                                 \
    spexit(stream);                                                \
  }

// and of two conditions for one variable
#define check_econd21(a, sign1_b1_sign0, sign2_b2, stream)              \
  if (a sign1_b1_sign0 a sign2_b2) {                                    \
    stream << "ERROR:\n";                                               \
    stream << '\n' << #a << #sign1_b1_sign0 << #a << #sign2_b2 << '\n'; \
    stream << #a << '=' << (a) << '\n';                                 \
    spexit(stream);                                                     \
  }

// and of two conditions for one variable
#define check_econd23(a, sign1, b1, sign0, sign2, b2, stream)               \
  if (a sign1 b1 sign0 a sign2 b2) {                                        \
    stream << "ERROR:\n";                                                   \
    stream << '\n'                                                          \
           << #a << #sign1 << #b1 << #sign0 << #a << #sign2 << #b2 << '\n'; \
    stream << #a << '=' << (a) << ' ' << #b1 << '=' << (b1) << ' ' << #b2   \
           << '=' << (b2) << '\n';                                          \
    spexit(stream);                                                         \
  }

// two conditions for four variables
#define check_econd24(a1, sign1, b1, sign0, a2, sign2, b2, stream)            \
  if (a1 sign1 b1 sign0 a2 sign2 b2) {                                        \
    stream << "ERROR:\n";                                                     \
    stream << '\n'                                                            \
           << #a1 << #sign1 << #b1 << #sign0 << #a2 << #sign2 << #b2 << '\n'; \
    stream << #a1 << '=' << (a1) << ' ' << #b1 << '=' << (b1) << '\n';        \
    stream << #a2 << '=' << (a2) << ' ' << #b2 << '=' << (b2) << '\n';        \
    spexit(stream);                                                           \
  }

namespace Heed {

class ExcFromSpexit {
 public:
  ExcFromSpexit() {}
};

void spexit_action(std::ostream& file);
extern int s_throw_exception_in_spexit;  // if == 1, does exit(1) of abort()
                                         // depending on the key below
extern int s_exit_without_core;          // the key above have larger priority

}  // namespace Heed

// Normal exit:
#define spexit(stream)                                                   \
  {                                                                      \
    stream << "File is " << __FILE__ << " , line number is " << __LINE__ \
           << '\n';                                                      \
    spexit_action(stream);                                               \
  }

#endif
