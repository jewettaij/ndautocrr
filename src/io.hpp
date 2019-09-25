///   @file io.hpp
///   @brief  A series of (arguably unnecessary) functions that 
///           read strings and numbers from a file.
///   @date 2007-4-13

#ifndef _IO_HPP
#define _IO_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
using namespace std;
#include "err.hpp"



/// When a syntax error occurs in a file, which file did it occur in?
string     g_filename;

/// where in the file (did the error occur):
long long  g_line;    

/// comment indicator
char        g_comments_begin_with ='#';

/// comment terminator
char        g_comments_end_with = '\n';


/// @brief  Does character c belong to aC?
template<class C>
inline bool
BelongsToCstring(C c, C const *aC, C terminating_char)
{
  assert(aC);
  while(*aC != terminating_char)
  {
    if (c == *aC) return true;
    ++aC;
  }
  return false;
}


/// @brief  Skip over some characters if they are present at this location
///         in the file.

inline void
Skip(istream &in, 
     char const *aSkipTheseChars)
{
  assert(aSkipTheseChars);
  if (! in) {
    stringstream err_msg;
    err_msg << "Error in input: \"" << g_filename <<"\"\n"
      "      near line " << g_line <<": File ends prematurely." << endl;
    throw InputErr(err_msg.str().c_str());
  }

  char c;
  while(in.get(c))
  {
    if (c=='\n') g_line++;  // <-- keep track of what line we are on

    else if (c == g_comments_begin_with)  
    {
      // skip past this comment (do not assume single-line comments)
      do {
        in.get(c);
        if (in && (c=='\n'))  g_line++;
      } while (in && (c != '\0') && (c != EOF) && (c!=g_comments_end_with));
    }
    // Note: this code does not work with nested comments. fix this later..

    if (in && (! BelongsToCstring(c, aSkipTheseChars, '\0')))
    {
      in.putback(c);
      if (c=='\n') g_line--;
      break;
    }
  }
} // Skip()



/// @brief  Read a token from the stream.
///         (Tokens are words delimited by one of the characters
///          in the "terminators" argument, which are typically whitespace.)

inline void
ReadString(istream &in, 
           string &dest,
           char const *terminators)
{
  assert(terminators);
  if (! in) {
    stringstream err_msg;
    err_msg << "Error in input: \"" << g_filename <<"\"\n"
      "      near line " << g_line <<": File ends prematurely." << endl;
    throw InputErr(err_msg.str().c_str());
  }

  char c;
  while(in.get(c))
  {
    if (BelongsToCstring(c, terminators, '\0') ||
        (c == g_comments_begin_with))
    {
      in.putback(c);
      break;
    }
    else 
    {
      dest.push_back(c);
      if (c=='\n') 
        g_line++; //keep track of the number of newlines read so far.
    }
  }
} // ReadString()



/// @brief  Read a number from this location in the file.

template<class Scalar>
bool
ReadScalar(istream &in, 
           Scalar& dest,
           char const *ac_terminators,
           string *ps_dest,
           string::const_iterator *pstopped_at=nullptr)
{
  string s;
  ReadString(in, s, ac_terminators);

  // If requested, make a copy this string and return it to the caller
  if (ps_dest != nullptr)
    *ps_dest=s;
  if (pstopped_at != nullptr)
    *pstopped_at = ps_dest->begin();

  if (s.size() != 0)
  {
    // I would prefer to use the standard ANSI C function: strtod()
    // to parse the string and convert it to a floating point variable.
    // But I have to copy this text into a c_string to get around strtod()'s
    // syntax.  strtod() requires an argument of type "char *".  I could
    // use s.c_str(), but this is a pointer of type "const char *".
    // So, I need to copy the contents of s.c_str() into a temporary array,
    // and then invoke strtod() on it.
    char *ac = new char [s.size() + 1];
    strcpy(ac, s.c_str());
    char *pstart = ac;
    char *pstop;

    #ifdef STRTOLD_UNSUPPORTED
    dest = strtod(pstart, &pstop);
    #else
    dest = strtold(pstart, &pstop);//Useful but not standard ANSI C
    #endif
    // Now pstop points past the last valid char

    // If requested, inform the caller where the parsing stopped
    if (pstopped_at != nullptr)
      *pstopped_at = ps_dest->begin() + (pstop - pstart);

    delete [] ac;
    return (pstop - pstart == s.size()); //did parsing terminate prematurely?
  }
  else
    return false; //no number was successfully read

} //ReadScalar()



/// @brief  Read a number from this location in the file.
/// @overloaded

template<class Scalar>
Scalar
ReadScalar(istream &in, 
           char const *ac_terminators)
{
  Scalar dest;
  string s;
  if (! ReadScalar(in, dest, ac_terminators, &s))
  {
    stringstream err_msg;
    err_msg << 
      "Error in input: \"" << g_filename << "\"\n"
      "      near line " << g_line;
    if (! s.empty())
      cerr << ": \""<<s<<"\"\n";
    err_msg << 
      "      Expected a number." << endl;
    throw InputErr(err_msg.str().c_str());
  }
  return dest;
} //ReadScalar()




/// @brief  Read an integer from this location in the file.

bool
ReadInt(istream &in, 
       long long& dest,
       char const *ac_terminators,
       string *ps_dest,
       string::const_iterator *pstopped_at = nullptr)
{
  string s;
  ReadString(in, s, ac_terminators);

  // If requested, make a copy this string and return it to the caller
  if (ps_dest != nullptr)
    *ps_dest=s;
  if (pstopped_at != nullptr)
    *pstopped_at = ps_dest->begin();

  if (s != "")
  {
    // I would prefer to use the standard ANSI C function: strtod()
    // to parse the string and convert it to a floating point variable.
    // But I have to copy this text into a c_string to get around strtod()'s
    // syntax.  strtod() requires an argument of type "char *".  I could
    // use s.c_str(), but this is a pointer of type "const char *".
    // So, I need to copy the contents of s.c_str() into a temporary array,
    // and then invoke strtod() on it.
    char *ac = new char [s.size() + 1];
    strcpy(ac, s.c_str());
    char *pstart = ac;
    char *pstop;

    dest = strtol(pstart, &pstop, 10);
    // Now pstop points past the last valid char

    // If requested, inform the caller where the parsing stopped
    if (pstopped_at != nullptr)
      *pstopped_at = ps_dest->begin() + (pstop - pstart);

    delete [] ac;
    return (pstop - pstart == s.size()); //did parsing terminate prematurely?
  }
  else
    return false; //no number was successfully read

} //ReadInt()


/// @brief  Read an integer from this location in the file.
/// @overloaded

long long
ReadInt(istream &in, 
        char const *ac_terminators)
{
  long long dest;
  string s;
  if (! ReadInt(in, dest, ac_terminators, &s))
  {
    stringstream err_msg;
    err_msg << 
      "Error in input: \"" << g_filename << "\"\n"
      "      near line " << g_line;
    if (! s.empty())
      cerr << ": \""<<s<<"\"\n";
    err_msg << 
      "      Expected an integer." << endl;
    throw InputErr(err_msg.str().c_str());
  }
  return dest;
} //ReadInt()



#endif //#ifndef _IO_HPP
