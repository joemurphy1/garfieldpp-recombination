#ifndef PRSTREAM_H
#define PRSTREAM_H
/*
This is the main file which determines the output matters:
default streams and indentation.

There are two default streams in C++: cout and cerr.
In the program we often need to use stream for regular output and
the same or another stream for emergency or exstraordinary cases,
exceptions, error, etc. These logical streams can be realized not only to
default tty, but to files. Perhaps there are many ways to control this,
but the simplest one is the use symbolic stream
notations mcout and mcerr (my cout and my cerr) throughout the program,
and to bound them with real streams through trivial macro-driven replacements,
as done below.

The practice shows that whatever advanced debugger and proficient skills
the programmer has in his computer, sooner or later
he will encode the print of all members of each significant
not trivial class of his program in readable and understandable form.
At least it is so in numerical calculations.
Each such printing is usually controlled by a key determining
the level of details. At large level the user expects to see
the output of all structured components of the current object.
Then the initial call of object->print(stream, key)
triggers similar calls of print of components, usually with less
key of details, component->print(stream, key-1).
It is very important to allow the reader of such listing
to distinguish visually the output from main object,
from components of the main object, components of components, etc.
Also it is useful to emphasis the name of classes,
the printing of sequences of similar elements
such as elements of arrays performed in loops. The user may want
to stress any other structures appearing in his classes.
This is possible by inclusion of indentation with which
the objects are printed. There should be some number of blanks
established which should be printed prior to content of
any line printed from any object. Any object should be allowed to
add additional blanks and required to remove the additions at the
end of its output.

It appeared that it is not trivial to arrange such system that makes this
and is completely safe from any misuse. Despite of all the power of C++
it appears to be not possible without significant intrusion in internal
functioning of streams. There was some discussion in a news-group which
does not point to any appropriate solution. Therefore here this is done
by means which could be crititized in some respects by the lovers of
object-oriented approach, but it has the pronounced advantages that
it is compatible with any streams, it is convenient enough for practice,
it really works, and it really exists.

There is a class indentation and the global object of this
class called shortly "indn". This object keeps the current number of blanks
needed to insert in output listing  before each line.
This object is not tied to certain stream. Therefore this current
number will be valid for any stream.

The indentation is invoked if you print
  Imcout<<something     // indentation is invoked
instead of
  mcout<<something      // indentation is not invoked
Also you can use Ifile instead of file.
Also indentation in 2 sequencial lines will be made by
  Imcout<<something_in_1_line<<'\n'<<indn<<something_in_next_line<<'\n';
To change the number of blanks you put
  indn.n+=2;  // or 1, or 4, or what you want.
And don't forget to restore the previous value by
  indn.n-=2;
If you output the composite class with redefined operator<<
and this latter operator uses indentation, particularly it
starts from Ifile<<...,
and you want this class continues the line, for example:
  Ifile<<"name_of_object="<<object;
then you probably don't want to see additional blanks between
"name_of_object=" and the first line of the object itself.
Then you want to prohibit the indentation when printing the first line of
the object. To provide this you can print:
  Ifile<<"name_of_object="<<noindent<<object;
The directive "noindent" turns off the request for indentation
immidiately after it, but switched off later. Thus indentation will be
automatically turned on after when the first file<<indn is met and skipped.
The rest of the object will be printed with correct indentation despite of
"noindent" statement.

In addition, at the bottom of this file a few useful Iprint-like macro are
defined. The idea is to print not only variable, but start from its name
and "=". Thus instead of
  Imcout << "my_favorite_variable=" << my_favorite_variable << '\n';
you can use the shorter
  Iprintn(mcout, my_favorite_variable);

  Iprint - just print a variable without saying << '\n' at the end.
  Iprintn - print a variable and pass to next line.

Copyright (c) 2001 I. B. Smirnov

Permission to use, copy, modify, distribute and sell this file
and its documentation for any purpose is hereby granted without fee,
provided that the above copyright notice, this permission notice,
and notices about any modifications of the original text
appear in all copies and in supporting documentation.
It is provided "as is" without express or implied warranty.
*/

#include <iostream>

namespace Heed {

#define Iprint(file, name) file << #name << "=" << name;
#define Iprintn(file, name) file << #name << "=" << name << '\n';
// addition is convenient as notation of units

#define Iprint2n(file, name1, name2) \
  file << #name1 << "=" << name1 << ", " << #name2 << "=" << name2 << '\n';
#define Iprint3n(file, name1, name2, name3)                                \
  file << #name1 << "=" << name1 << ", " << #name2 << "=" << name2 << ", " \
       << #name3 << "=" << name3 << '\n';
#define Iprint4n(file, name1, name2, name3, name4)                         \
  file << #name1 << "=" << name1 << ", " << #name2 << "=" << name2 << ", " \
       << #name3 << "=" << name3 << ", " << #name4 << "=" << name4 << '\n';

// simultaneously for all classes. Useful for writing "persistence classes"
// by standard <</>> operators.
// If instead of this one tries to use special functions like
// class_name::short_write,
// he finds an obstacle that such functions are absent for inbuilt types.
}  // namespace Heed

#endif
