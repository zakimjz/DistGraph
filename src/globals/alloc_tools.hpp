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

#ifndef __ALLOC_TOOLS_HPP__
#define __ALLOC_TOOLS_HPP__

#include <memory_checker.hpp>
#include <logger.hpp>
#include <vector>
#include <cassert>
#include <new>


#define NEW_INT_ARRAY(__ptr__, __size__, __logger__) {                 \
    try{ \
      __ptr__ = new int[__size__];  \
    }catch(std::bad_alloc& ba) { \
      CRITICAL_ERROR(__logger__, "bad alloc error at " << __FILE__ << "@" << __LINE__ << "; err string: " << ba.what()); \
      throw std::runtime_error("Error in memory allocation"); \
    } \
    REGISTER_MBLOCK2(__ptr__, __size__);          \
}

#define DELETE_INT_ARRAY(__ptr__, __logger__) {             \
    try{ \
      delete[] __ptr__;  \
    }catch(...) { \
      CRITICAL_ERROR(__logger__, "dealloc error at " << __FILE__ << "@" << __LINE__ ); \
      throw std::runtime_error("Error in memory deallocation"); \
    } \
    UNREGISTER_MBLOCK(__ptr__); \
}

#endif

