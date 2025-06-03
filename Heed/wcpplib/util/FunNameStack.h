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
#ifdef FUNNAMESTACK
#define mfunname(string)                   \
  static const char* FunNameIIII = string; \
  FunNameWatch funnw(FunNameIIII)
#else
#define mfunname(string)
#endif

// Permanent definitions
#define mfunnamep(string)                  \
  static const char* FunNameIIII = string; \
  FunNameWatch funnw(FunNameIIII)

// Switch on/off checks
#define DO_CHECKS
#ifdef DO_CHECKS

#ifdef FUNNAMESTACK

#define check_econd1(condition, a1, stream) \
  if (condition) {                          \
    funnw.ehdr(stream);                     \
    stream << '\n' << #condition << '\n';   \
    stream << #a1 << '=' << (a1) << '\n';   \
    spexit(stream);                         \
  }
#define check_wcond1(condition, a1, stream) \
  if (condition) {                          \
    funnw.whdr(stream);                     \
    stream << '\n' << #condition << '\n';   \
    stream << #a1 << '=' << (a1) << '\n';   \
  }

#define check_econd11(a, signb, stream)     \
  if (a signb) {                            \
    funnw.ehdr(stream);                     \
    stream << '\n' << #a << #signb << '\n'; \
    stream << #a << '=' << (a) << '\n';     \
    spexit(stream);                         \
  }

#define check_econd12(a, sign, b, stream)                          \
  if (a sign b) {                                                  \
    funnw.ehdr(stream);                                            \
    stream << '\n' << #a << #sign << #b << '\n';                   \
    stream << #a << '=' << (a) << ' ' << #b << '=' << (b) << '\n'; \
    spexit(stream);                                                \
  }

// condition + additional any commands
#define check_econd11a(a, signb, add, stream) \
  if (a signb) {                              \
    funnw.ehdr(stream);                       \
    stream << '\n' << #a << #signb << '\n';   \
    stream << #a << '=' << (a) << '\n';       \
    stream << add;                            \
    spexit(stream);                           \
  }

#define check_econd12a(a, sign, b, add, stream)                    \
  if (a sign b) {                                                  \
    funnw.ehdr(stream);                                            \
    stream << '\n' << #a << #sign << #b << '\n';                   \
    stream << #a << '=' << (a) << ' ' << #b << '=' << (b) << '\n'; \
    stream << add;                                                 \
    spexit(stream);                                                \
  }

// and of two conditions for one variable
#define check_econd21(a, sign1_b1_sign0, sign2_b2, stream)              \
  if (a sign1_b1_sign0 a sign2_b2) {                                    \
    funnw.ehdr(stream);                                                 \
    stream << '\n' << #a << #sign1_b1_sign0 << #a << #sign2_b2 << '\n'; \
    stream << #a << '=' << (a) << '\n';                                 \
    spexit(stream);                                                     \
  }

// and of two conditions for one variable
#define check_econd23(a, sign1, b1, sign0, sign2, b2, stream)               \
  if (a sign1 b1 sign0 a sign2 b2) {                                        \
    funnw.ehdr(stream);                                                     \
    stream << '\n'                                                          \
           << #a << #sign1 << #b1 << #sign0 << #a << #sign2 << #b2 << '\n'; \
    stream << #a << '=' << (a) << ' ' << #b1 << '=' << (b1) << ' ' << #b2   \
           << '=' << (b2) << '\n';                                          \
    spexit(stream);                                                         \
  }

// two conditions for four variables
#define check_econd24(a1, sign1, b1, sign0, a2, sign2, b2, stream)            \
  if (a1 sign1 b1 sign0 a2 sign2 b2) {                                        \
    funnw.ehdr(stream);                                                       \
    stream << '\n'                                                            \
           << #a1 << #sign1 << #b1 << #sign0 << #a2 << #sign2 << #b2 << '\n'; \
    stream << #a1 << '=' << (a1) << ' ' << #b1 << '=' << (b1) << '\n';        \
    stream << #a2 << '=' << (a2) << ' ' << #b2 << '=' << (b2) << '\n';        \
    spexit(stream);                                                           \
  }

#else  // without FUNNAMESTACK, print only condition

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

#endif

#else  // without checks

#define check_econd1(condition, a1, stream)
#define check_wcond1(condition, a1, stream)

#define check_econd11(a, signb, stream)
#define check_econd12(a, sign, b, stream)
#define check_econd11a(a, signb, add, stream)
#define check_econd12a(a, sign, b, add, stream)
// and of two conditions for one variable
#define check_econd21(a, sign1_b1_sign0, sign2_b2, stream)
// and of two conditions for one variable
#define check_econd23(a, sign1, b1, sign0, sign2, b2, stream)
// two conditions for four variables
#define check_econd24(a1, sign1, b1, sign0, a2, sign2, b2, stream)

#endif

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
    stream << FunNameStack::instance();                                  \
    stream << "File is " << __FILE__ << " , line number is " << __LINE__ \
           << '\n';                                                      \
    spexit_action(stream);                                               \
  }

namespace Heed {

class FunNameStack {
 public:
  static FunNameStack& instance();

  // make instances of the class uncopyable to hidden exsemplar,
  // but copyable away.
  //

  FunNameStack(const FunNameStack& f);
  FunNameStack& operator=(const FunNameStack& f);

  FunNameStack(void);  // usually called inly from instance()
 private:
  int qname;
  static const int pqname = 1000;
  char* name[pqname];
  int s_init;  // 1 - sign that it is inited,
               // any other value - not inited
               // Now this variable is used only for debug printing.

  int s_act;  // 1 - active, 0 - not active
              //( to switch without recompilation)
 public:
  int s_print;  // works only is s_act=1;
                // 0 - no print
                // 1 - at each entry print the name of current function
                // 2 - at each entry and each exit
                //    print the name of current function
                // 3 - at each entry print all stack
                // 4 - at each entry and exit print all stack
 private:
  int nmode;  // 0 - name points to a string in outside world
              //     used for global object funnamestack
              // 1 - name is inited by new and deleted by delete
              //     may be used for sending as parameter of exception
              //     It is not convenient for normal work due
              //     to time expense

  std::ostream& printname(std::ostream& file, int n);

 public:
  void set_parameters(int fs_act = 1, int fs_print = 0);
  ~FunNameStack();

 private:
  void printput(std::ostream& file);  // called at insertion of new name in
                                      // stack
  void printdel(std::ostream& file);  // called at deletion of name from stack
 public:
  inline int put(const char* fname) {
    // if(s_init != 1) init();
    if (s_act != 1) return 0;
    if (qname >= pqname) {
      mcerr << "FunNameStack::put: error: qname == pqname\n";
      mcerr << "Most oftenly this happens due to infinite recursion.\n";
      mcerr << "*this=" << (*this);
      exit(1);
    }
    name[qname++] = const_cast<char*>(fname);
    if (s_print > 0) {
      printput(mcout);
    }
    return qname - 1;
  }
  inline void del(int nname) {
    if (s_act != 1) return;
    if (nname != qname - 1) {
      // not last
      qname = nname;
    } else {
      if (s_print > 0) {
        printdel(mcout);
      }
      qname--;
    }
  }
  inline void replace(const char* fname) {
    // if(s_init != 1) init();
    if (s_act != 1) return;
    if (qname >= pqname) {
      mcerr << "FunNameStack::put: error: qname == pqname\n";
      mcerr << "Most oftenly this happens due to infinite recursion.\n";
      mcerr << "*this=" << (*this);
      exit(1);
    }
    name[qname - 1] = const_cast<char*>(fname);
    if (s_print > 0) printput(mcout);
  }
  friend std::ostream& operator<<(std::ostream& file, const FunNameStack& f);
};
std::ostream& operator<<(std::ostream& file, const FunNameStack& f);

// Object of this class FunNameWatch is created at the beginning of function.
// At creation it writes function name in stack.
// At deletion it deletes the name from stack,
// checking that it is last.

class FunNameWatch {
  int nname;         // index of this name in array of FunNameStack
  const char* name;  // it is memorized independenlty on s_act.
                     // Used for printing of headers.
 public:
  inline FunNameWatch(const char* fname) {
    nname = FunNameStack::instance().put(fname);
  }
  inline ~FunNameWatch() {
    if (nname >= 0) FunNameStack::instance().del(nname);
  }

  // print header
  std::ostream& hdr(std::ostream& file) const {
    file << name << ": ";
    return file;
  }
  // print header with word WARNING
  std::ostream& whdr(std::ostream& file) const {
    file << name << ": WARNING:\n";
    return file;
  }
  // print header with word ERROR
  std::ostream& ehdr(std::ostream& file) const {
    file << name << ": ERROR:\n";
    return file;
  }
};
std::ostream& operator<<(std::ostream& file, const FunNameWatch& f);

}  // namespace Heed

#endif
