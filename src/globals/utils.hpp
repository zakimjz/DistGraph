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

#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <string>
#include <sys/time.h>
#include <time.h>
#include <vector>
#include <map>

#define MEGABYTES (1024 * 1024)

namespace utils {

std::string trim(const std::string& s, const std::string& drop = " ");


bool is_na_symbol(const char *str);
bool is_symbol(const char *str);
bool is_string(const char *cstr);

bool parse_int(const char *str, int &ret);
bool parse_uint(const char *cstr, unsigned int &ret);
bool parse_double(const char *str, double &ret);
bool parse_long(const char *str, long &ret);
bool parse_ulong(const char *str, unsigned long &ret);

bool is_int(const char *str);
bool is_double(const char *str);
bool is_uint(const char *cstr);
bool is_long(const char *cstr);
bool is_ulong(const char *cstr);


bool is_symbol(const char *str);
bool is_na_symbol(const char *str);

void split(const std::string& str,
           std::vector<std::string>& tokens,
           const std::string& delimiters = " ");


time_t get_sec(timeval & start, timeval &end);
time_t get_usec(timeval & start, timeval &end);
double get_time_diff(timeval & start, timeval &end);

long get_program_size();
long get_rss();
long get_rss_limit();

std::string get_current_timestamp();
std::string print_vector(const std::vector<int> &vec);
std::string print_array(int *a, int len);
std::string print_map(const std::map<int, int> &imap);

unsigned long hash(std::string str);

} // namespace utils
#endif
