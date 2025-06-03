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

}  // namespace Heed
