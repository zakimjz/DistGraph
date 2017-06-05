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

#ifndef __MEMORY_CHECKER_HPP__
#define __MEMORY_CHECKER_HPP__

#include <string>
#include <sstream>
#include <stdint.h>

#define REGISTER_MBLOCK(__ptr__)    memory_checker::register_memory_block(__ptr__, __LINE__, __FILE__)
#define UNREGISTER_MBLOCK(__ptr__)  memory_checker::unregister_memory_block(__ptr__)

#define REGISTER_MBLOCK2(__ptr__, __size__)    memory_checker::register_memory_block(__ptr__, __size__, __LINE__, __FILE__)

namespace memory_checker {

struct allocated_block_descriptor {
  allocated_block_descriptor(int *&ptr, int line, std::string file) {
    this->line = line;
    this->file = file;
    this->ptr = ptr;
    size = 0;
  }
  allocated_block_descriptor(int *&ptr, unsigned int size, int line, std::string file) {
    this->line = line;
    this->file = file;
    this->ptr = ptr;
    this->size = size;
  }
  int line;
  std::string file;
  int *ptr;
  unsigned int size;

  std::string to_string() const {
    std::stringstream ss;
    ss << file << "@" << line << "; ptr: " << ptr << "; size: " << size;
    return ss.str();
  }
};
void register_memory_block(int *&ptr, int line, const std::string &file);
void register_memory_block(int *&ptr, unsigned int size, int line, const std::string &file);
void unregister_memory_block(int *&ptr);
void print_unallocated_blocks(bool report_leaks = true);
void detect_memory_leaks(bool report_leaks = true);
bool is_block_allocated(int *&ptr);
std::string print_memory_usage();
double get_memory_usage();
double get_max_memory_usage_mb();

} // namespace memory_checker

#endif


