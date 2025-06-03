/*
Copyright (c) 1999-2003 I. B. Smirnov

Permission to use, copy, modify, distribute and sell this file for any purpose
is hereby granted without fee, provided that the above copyright notice,
this permission notice, and notices about any modifications of the original
text appear in all copies and in supporting documentation.
The file is provided "as is" without express or implied warranty.
*/

#include "wcpplib/util/FunNameStack.h"

#include <cstring>
#include <iomanip>
#include <iostream>

namespace Heed {

int s_throw_exception_in_spexit = 0;
int s_exit_without_core = 0;

FunNameStack& FunNameStack::instance() {
  // static FunNameStack inst;  // According to some Internet site this is
  // Meyer's approach,
  // // time of descruction is not determined.
  // // So it can potentially be destructed before it is still in use.
  // return inst;

  static FunNameStack* inst = NULL;  // According to this site it is
                                     // "GoF" approach,
                                     // destruction is not performed at all.
#ifdef USE_TOGETHER_WITH_CLEAN_NEW
#if defined(MAINTAIN_KEYNUMBER_LIST) && defined(USE_BOOST_MULTITHREADING)
  MemoriseIgnore::instance().ignore();
#else
  s_ignore_keynumberlist = 1;  // to avoid report from delete at deletion
#endif
#endif
  if (!inst) inst = new FunNameStack();
#ifdef USE_TOGETHER_WITH_CLEAN_NEW
#if defined(MAINTAIN_KEYNUMBER_LIST) && defined(USE_BOOST_MULTITHREADING)
  MemoriseIgnore::instance().not_ignore();
#else
  s_ignore_keynumberlist = 0;
#endif
#endif
  return *inst;
}

FunNameStack::FunNameStack(void)
    :  // usually called inly from instance()
      s_init(1),
      s_act(1),
      s_print(0),
      nmode(0) {
  qname = 0;
  for (int n = 0; n < pqname; n++) name[n] = NULL;
}

std::ostream& FunNameStack::printname(std::ostream& file, int n)  //
{
  file << name[n];
  return file;
}

void spexit_action(std::ostream& file) {
  file << "spexit_action: the streams will be now flushed\n";
  file.flush();
  mcout.flush();
  mcerr.flush();
  if (s_throw_exception_in_spexit != 1) {
    if (s_exit_without_core == 1) {
      file << "spexit_action: the exit(1) function is called\n";
      exit(1);
    } else {
      file << "spexit_action: the abort function is called\n";
      abort();
    }
  } else {
    file << "spexit_action: an exception is now called\n";
    throw ExcFromSpexit();
  }
}

FunNameStack::FunNameStack(const FunNameStack& f) { *this = f; }

FunNameStack& FunNameStack::operator=(const FunNameStack& f) {
  if (this == &(FunNameStack::instance())) {
    mcerr << "ERROR in FunNameStack& FunNameStack::operator=(const "
             "FunNameStack& f)\n";
    mcerr << "Attempt to copy to operative FunNameStack\n";
    mcerr << "You don't need and should not initialize main FunNameStack "
             "directly\n";
    mcerr << "If you want to change its parameters, use "
             "FunNameStack::instance()\n";
    Iprintn(mcout, this);
    Iprintn(mcout, &(FunNameStack::instance()));
    Iprintn(mcout, (*this));
    Iprintn(mcout, (FunNameStack::instance()));

    spexit(mcerr);
  }
  qname = f.qname;
  s_init = f.s_init;
  s_act = f.s_act;
  s_print = f.s_print;
  nmode = f.nmode;
  for (int n = 0; n < f.qname; n++) {
    if (nmode == 0) {
      name[n] = f.name[n];
    } else {
      int l = strlen(f.name[n]) + 1;
#ifdef USE_TOGETHER_WITH_CLEAN_NEW
      s_ignore_keynumberlist = 1;
#endif
      name[n] = new char[l];
      strcpy(name[n], f.name[n]);
#ifdef USE_TOGETHER_WITH_CLEAN_NEW
      s_ignore_keynumberlist = 0;
#endif
    }
  }
  return *this;
}

void FunNameStack::set_parameters(int fs_act, int fs_print) {
  s_act = fs_act;
  s_print = fs_print;  // only to correct this
}

FunNameStack::~FunNameStack() {
#ifdef USE_TOGETHER_WITH_CLEAN_NEW
  s_ignore_keynumberlist = 1;
#endif
  if (nmode == 1)
    for (int n = 0; n < qname; n++) delete name[n];
#ifdef USE_TOGETHER_WITH_CLEAN_NEW
  s_ignore_keynumberlist = 0;
#endif
}

void FunNameStack::printput(std::ostream& file) {
  if (s_print == 1 || s_print == 2) {
    file << "FunNameStack::put: qname =" << qname << " last name=";
    printname(file, qname - 1) << '\n';
  } else if (s_print >= 3) {
    file << "FunNameStack::put:\n" << (*this);
  }
}
void FunNameStack::printdel(std::ostream& file) {
  if (s_print == 2) {
    file << "FunNameStack::del: qname =" << qname << " last name=";
    printname(file, qname - 1) << '\n';
  } else if (s_print == 4) {
    file << "FunNameStack::del:\n" << (*this);
  }
}

std::ostream& operator<<(std::ostream& file, const FunNameStack& f) {
  if (f.s_act == 1) {
    file << "FunNameStack: s_init=" << f.s_init << " qname=" << f.qname << '\n';
    for (int n = 0; n < f.qname; n++) {
      file << std::setw(3) << n << "  " << f.name[n] << " \n";
    }
  }
  return file;
}

std::ostream& operator<<(std::ostream& file, const FunNameWatch& f) {
  f.hdr(file);
  return file;
}

}  // namespace Heed
