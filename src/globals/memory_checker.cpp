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

#include <memory_checker.hpp>
#include <map>
#include <logger.hpp>
#include <stdexcept>

namespace memory_checker {

typedef std::map<uintptr_t, allocated_block_descriptor> allocated_memory_blocks_t;
static allocated_memory_blocks_t mem_blocks;

static unsigned int current_allocated_size = 0;
static unsigned int max_allocated_size = 0;

void register_memory_block(int *&ptr, unsigned int size, int line, const std::string &file)
{
  if(ptr == 0) {
    return;
  }
  //std::cout<<"registering memory block with ptr: " << ptr <<" size " << size << " from " << file << "@" << line << std::endl;
  if(mem_blocks.find((uintptr_t) ptr) != mem_blocks.end()) {
    CRITICAL_ERROR(*Logger::get_logger("MEMCHECK"),  "registering already registered memory block with ptr: " << ptr << " from " << file << "@" << line);
    abort();
  } // if
  mem_blocks.insert(std::make_pair((uintptr_t) ptr, allocated_block_descriptor(ptr, size, line, file)));
  current_allocated_size += size;
  if(max_allocated_size < current_allocated_size) max_allocated_size = current_allocated_size;

}

void register_memory_block(int *&ptr, int line, const string &file)
{
  if(ptr == 0) {
    return;
  }
  if(mem_blocks.find((uintptr_t) ptr) != mem_blocks.end()) {
    CRITICAL_ERROR(*Logger::get_logger("MEMCHECK"), "registering already registered memory block with ptr: " << ptr << " from " << file << "@" << line);
    CRITICAL_ERROR(*Logger::get_logger("MEMCHECK"), "already registered from: " << mem_blocks.find((uintptr_t) ptr)->second.to_string());
    abort();
  } // if
  mem_blocks.insert(std::make_pair((uintptr_t) ptr, allocated_block_descriptor(ptr, line, file)));
} // register_memory_block


void unregister_memory_block(int *&ptr)
{
  if(ptr == 0) {
    return;
  }
  //std::cout<<"unregistering memory block with ptr: " << ptr<<std::endl;
  std::map<uintptr_t, allocated_block_descriptor>::iterator it = mem_blocks.find((uintptr_t) ptr);
  if(it == mem_blocks.end()) {
    CRITICAL_ERROR(*Logger::get_logger("MEMCHECK"),  "un-registering memory block that was not registered with ptr: " << ptr);
    abort();
  } // if
  current_allocated_size -= it->second.size;
  mem_blocks.erase((uintptr_t) ptr);
}

std::ostream &operator<<(ostream &out, const allocated_block_descriptor &block) {
  out << "memory block, ptr: " << block.ptr << "; allocated at " << block.file << "@" << block.line << " of size: " << block.size;
  return out;
}


void print_unallocated_blocks(bool report_leaks)
{
  for(allocated_memory_blocks_t::iterator it = mem_blocks.begin(); it != mem_blocks.end(); it++) {
    cout << it->second << endl;
  } // for it
} // print_unallocated_blocks


void detect_memory_leaks(bool report_leaks)
{
  if(mem_blocks.empty()) {
    INFO(*Logger::get_logger("MEMCHECK"), "all memory blocks were freed !");
    return;
  }

  for(allocated_memory_blocks_t::iterator it = mem_blocks.begin(); it != mem_blocks.end(); it++) {
    cout << it->second << endl;
  } // for it
  throw std::runtime_error("Detected memory leaks in memory !!!");
} // print_unallocated_blocks


bool is_block_allocated(int *&ptr)
{
  return mem_blocks.find((uintptr_t) ptr) != mem_blocks.end();
}


double get_memory_usage()
{
  int total_size = 0;
  for(allocated_memory_blocks_t::iterator it = mem_blocks.begin(); it != mem_blocks.end(); it++) {
    //cout << it->second << endl;
    total_size += it->second.size;
  } // for it
  return double(total_size);
}

double get_max_memory_usage_mb()
{
  return double(max_allocated_size) / (1024.0L * 1024.0L);
}

std::string print_memory_usage()
{
  std::stringstream ss;
  ss << "total memory usage: " << double(get_memory_usage()) / (1024.0L * 1024.0L);
  return ss.str();
}


} // namespace memory_checker




