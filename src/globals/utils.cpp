/*
 * This source is part of the single graph mining algorithm.
 *
 * Copyright 2014 Robert Kessl
 * Copyright 2015-2016 Nilothpal Talukder
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <string>
#include <vector>
#include <algorithm>

#include <iostream>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <fstream>
#include <sys/resource.h>


#include <utils.hpp>
#include <logger.hpp>

using std::string;
using std::vector;


namespace utils {

std::string trim(const std::string& s, const std::string& drop)
{
  string r = s;
  r.erase(r.find_last_not_of(drop) + 1);
  return r.erase(0,r.find_first_not_of(drop));
}

/** 'str' must be trimmed !!!!
 */
bool is_double(const char *cstr)
{
  double d;
  return parse_double(cstr, d);
  /*
     char *ptr = 0;
     const char *cstr = str.c_str();
     uint_t cstr_len = strlen(cstr);
     strtod(cstr, ptr);
     // TODO: check space or trim str !!!
     if(ptr == cstr+cstr_len) return true;
     return false;
   */
}


bool parse_double(const char *cstr, double &ret)
{
  char *ptr = 0;
  //const char *cstr = str.c_str();
  unsigned int cstr_len = strlen(cstr);
  double d = strtod(cstr, &ptr);
  // TODO: check space or trim str !!!
  if(ptr == cstr + cstr_len) {
    ret = d;
    return true;
  }
  return false;
}




/* parse_XXX section
 */
bool parse_int(const char *cstr, int &ret)
{
  char *ptr = 0;
  //const char *cstr = str.c_str();
  unsigned int cstr_len = strlen(cstr);
  int i = (int) strtol(cstr, &ptr, 10);
  // TODO: check space or trim str !!!
  if(ptr == cstr + cstr_len) {
    ret = i;
    return true;
  }
  return false;
}



bool parse_uint(const char *cstr, unsigned int &ret)
{
  char *ptr = 0;
  //const char *cstr = str.c_str();
  unsigned int cstr_len = strlen(cstr);
  unsigned int i = (unsigned int) strtoul(cstr, &ptr, 10);
  // TODO: check space or trim str !!!
  if(ptr == cstr + cstr_len) {
    ret = i;
    return true;
  }
  return false;
}


bool parse_long(const char *cstr, long &ret)
{
  char *ptr = 0;
  unsigned int cstr_len = strlen(cstr);
  long i = strtol(cstr, &ptr, 10);
  if(ptr == cstr + cstr_len) {
    ret = i;
    return true;
  }
  return false;
}

bool parse_ulong(const char *cstr, unsigned long &ret)
{
  char *ptr = 0;
  unsigned int cstr_len = strlen(cstr);
  long i = strtol(cstr, &ptr, 10);
  if(ptr == cstr + cstr_len) {
    ret = i;
    return true;
  }
  return false;
}




/* is_XXX section
 */
bool is_int(const char *cstr)
{
  int i;
  return parse_int(cstr, i);
}


bool is_uint(const char *cstr)
{
  unsigned int i;
  return parse_uint(cstr, i);
}

bool is_long(const char *cstr)
{
  long i;
  return parse_long(cstr, i);
}


bool is_ulong(const char *cstr)
{
  unsigned long i;
  return parse_ulong(cstr, i);
}



/** str must be trimmed !!!
 */
bool is_symbol(const char *cstr)
{
  if(isdigit(cstr[0])) return false;
  if(isalpha(cstr[0])) {
    for(int i = 1; i < strlen(cstr); i++) {
      if(!isalnum(cstr[i]) && cstr[i] != '_') return false;
    }
    return true;
  }
  return false;
}


bool is_na_symbol(const char *str)
{
  if(strcmp(str, "NA") == 0) return true;
  if(strcmp(str, "na") == 0) return true;
  if(strcmp(str, "N/A") == 0) return true;
  if(strcmp(str, "n/a") == 0) return true;

  return false;
}


bool is_string(const char *cstr) {
  unsigned int l = strlen(cstr);
  for(int i = 0; i < l; i++) {
    if(isprint(cstr[i]) == 0) return false;
  }
  return true;
}


void split(const string& str,
           vector<string>& tokens,
           const string& delimiters)
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}












time_t get_sec(timeval & start, timeval &end)
{
  time_t sec;
  suseconds_t usec;

  sec = end.tv_sec - start.tv_sec;
  usec = end.tv_usec - start.tv_usec;
  if(usec < 0) {
    sec--;
    usec = 1000000 + usec;
  }

  return sec;
}

time_t get_usec(timeval & start, timeval &end)
{
  time_t sec;
  suseconds_t usec;

  sec = end.tv_sec - start.tv_sec;
  usec = end.tv_usec - start.tv_usec;
  if(usec < 0) {
    sec--;
    usec = 1000000 + usec;
  }

  return usec;
}


double get_time_diff(timeval & start, timeval &end)
{
  double time = (double) end.tv_sec + (double) end.tv_usec / 1000000.0L
                - ((double) start.tv_sec + (double) start.tv_usec / 1000000.0L);
  return time;
}

std::string get_current_timestamp()
{
  char strtime[100];
  tm *local_time = 0;
  time_t unix_time = time(0);

  local_time = localtime(&unix_time);
  strftime(strtime, 100, "%z", local_time);
  free(local_time);
  return string(strtime);
}


long get_program_size()
{
  long pg_size;
  pg_size = (long) getpagesize();
  std::ifstream proc;
  proc.open("/proc/self/statm", std::ios::in);
  long tmp;
  proc >> tmp;

  return tmp * pg_size;
}


long get_rss()
{
  long pg_size;
  pg_size = (long) getpagesize();
  std::ifstream proc;
  proc.open("/proc/self/statm", std::ios::in);
  long tmp;
  proc >> tmp;
  proc >> tmp;

  return tmp * pg_size;
}



long get_rss_limit()
{
  rlimit l;
  getrlimit(RLIMIT_AS, &l);
  long pg_size;
  pg_size = getpagesize();
  if(l.rlim_max == RLIM_INFINITY) return -1;
  return long(l.rlim_max) * pg_size;
}

std::string print_map(const std::map<int, int> &imap){
  std::stringstream ss;
  for(std::map<int, int>::const_iterator it = imap.begin(); it != imap.end(); it++) {
    ss << it->first << " : " << it->second << " , ";
  }
  return ss.str();
}

std::string print_vector(const std::vector<int> &vec) {
  if(vec.empty() == true) return "";
  std::stringstream ss;
  ss << vec[0];
  for(int i = 1; i < vec.size(); i++) ss << " " << vec[i];
  return ss.str();
}

std::string print_array(int* a, int len) {
  if(a == 0 or len == 0) return "";
  std::stringstream ss;
  ss << a[0];
  for(int i = 1; i < len; i++) ss << " " << a[i];
  return ss.str();
}

unsigned long hash(std::string str)
{
  unsigned long hash = 7919;
  int c;

  const char *cstr = str.c_str();
  while (c = *cstr++)
    hash = ((hash << 5) + hash) + c;     /* hash * 33 + c */

  return hash;
}

} // namespace utils
