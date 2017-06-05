/*
 * This source is part of the single graph mining algorithm.
 *
 * Copyright 2004 Taku Kudo
 * Copyright 2014 Robert Kessl
 * Copyright 2016 Nilothpal Talukder
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

#include <types.hpp>
#include <dfs_code.hpp>
#include <graph_types.hpp>
#include <sstream>
#include <stdexcept>
#include <cstdio>
#include <string>
#include <utils.hpp>
#include <exceptions.hpp>
#include <string.h>

using std::runtime_error;
using std::string;
using std::stringstream;

namespace types {


size_t DFS::get_serialized_size() const
{
  return sizeof(int) * 5;
}


size_t DFS::get_serialized_size(char *buffer, size_t buffer_size) const
{
  return *((int*)buffer);
}


size_t DFS::serialize(char *buffer, size_t buffer_size) const
{
  if(buffer_size < sizeof(int) * 5) throw runtime_error("Buffer too small");
  int *buf = (int*) buffer;
  buf[0] = from;
  buf[1] = to;
  buf[2] = fromlabel;
  buf[3] = elabel;
  buf[4] = tolabel;
  return 5 * sizeof(int);
}


size_t DFS::deserialize(char *buffer, size_t buffer_size)
{
  if(buffer_size < get_serialized_size()) throw runtime_error("Buffer too small");
  int *buf = (int*) buffer;
  from = buf[0];
  to = buf[1];
  fromlabel = buf[2];
  elabel = buf[3];
  tolabel = buf[4];
  return 5 * sizeof(int);
}

size_t DFS::serialize(std::ostream &) const
{
  throw method_unimplemented("DFS::serialize");
}


size_t DFS::deserialize(std::istream &)
{
  throw method_unimplemented("DFS::serialize");
}


/**
 * RMPath is in the opposite order then the DFS code, i.e., the
 * indexes into DFSCode go from higher numbers to lower numbers.
 *
 *
 *
 */
const RMPath &DFSCode::buildRMPath()
{
  rmpath.clear();

  int old_from = -1;

  for(int i = size() - 1; i >= 0; --i) {
    if((*this)[i].from < (*this)[i].to &&  // forward
       (rmpath.empty() || old_from == (*this)[i].to)) {
      rmpath.push_back(i);
      old_from = (*this)[i].from;
    }
  }

  return rmpath;
}


DFS DFS::parse_from_string(const char *str_dfs_code)
{
  size_t len = strlen(str_dfs_code);
  int from;
  int to;
  int fromlabel;
  int elabel;
  int tolabel;
  char F_B;

  int r = sscanf(str_dfs_code, "%c(%d %d %d %d %d)", &F_B, &from, &to, &fromlabel, &elabel, &tolabel);
  if(r != 6) throw runtime_error("error: could not read dfs code from string");

  return DFS(from, to, fromlabel, elabel, tolabel);
} // parse_from_string


string DFS::to_string(bool print_edge_type) const
{
  std::stringstream ss;
  if(print_edge_type) {
    if(is_forward()) ss << "F";
    else ss << "B";
  }
  ss << "(" << from << " " << to << " " << fromlabel << " " << elabel << " " << tolabel << ")";
  return ss.str();
}

bool DFS::is_forward() const
{
  return from < to;
}


bool DFS::is_backward() const
{
  return from > to;
}


bool DFS_less_then::bckwrd_bckwrd_less(const DFS &d1, const DFS &d2) const
{
  return (d1.to < d2.to) ||
         (d1.to == d2.to && d1.elabel < d2.elabel);
}


bool DFS_less_then::frwrd_bckwrd_less(const DFS &d1, const DFS &d2) const
{
  if(d1.is_backward() && d2.is_forward()) return true;
  return false;
}


bool DFS_less_then::frwrd_frwrd_less(const DFS &d1, const DFS &d2) const
{
  if(d1.from > d2.from) return true;
  if(d1.from < d2.from) return false;
  if(d1.from == d2.from) {
    //bool tmp = (d1.elabel < d2.elabel);
    if(d1.elabel < d2.elabel) return true;
    if(d1.elabel > d2.elabel) return false;

    return (d1.tolabel < d2.tolabel);
  } // if
  return false;
}


bool DFS_less_then::operator()(const DFS &d1, const DFS &d2) const
{
  if(d1.is_backward() && d2.is_backward()) {
    bool result = bckwrd_bckwrd_less(d1, d2);
    return result;
  }
  if(d1.is_forward() && d2.is_forward()) {
    bool result = frwrd_frwrd_less(d1, d2);
    return result;
  }

  bool result = frwrd_bckwrd_less(d1, d2);
  return result;
}


std::ostream &DFSCode::write(std::ostream &os) const
{
  if(size() == 0) return os;

  os << "(" << (*this)[0].fromlabel << ") " << (*this)[0].elabel << " (0f" << (*this)[0].tolabel << ");";

  for(unsigned int i = 1; i < size(); ++i) {
    if((*this)[i].from < (*this)[i].to) {
      os << " " << (*this)[i].elabel << " (" << (*this)[i].from << "f" << (*this)[i].tolabel << ");";
    } else {
      os << " " << (*this)[i].elabel << " (b" << (*this)[i].to << ");";
    }
  }

  return os;
}


std::string DFSCode::to_string(bool print_edge_type) const
{
  if(empty()) return "";
  std::stringstream ss;
  //write(ss);

  int i = 0;
  ss << (*this)[i].to_string(print_edge_type);
  i++;

  for(; i < size(); ++i) {
    ss << ";" << (*this)[i].to_string(print_edge_type);
  }

  return ss.str();
}


DFSCode DFSCode::read_from_str(const string &str)
{
  std::vector<std::string> vec_str;
  utils::split(str, vec_str, ";");

  DFSCode result;

  for(int i = 0; i < vec_str.size(); i++) {
    DFS d = DFS::parse_from_string(vec_str[i].c_str());
    result.push_back(d);
  } // for i

  return result;
}



std::ostream &operator<<(std::ostream &out, const DFSCode &code) {
  out << code.to_string();
  return out;
}


bool DFSCode::toGraph(Graph &g) const
{
  g.clear();

  for(DFSCode::const_iterator it = begin(); it != end(); ++it) {
    g.resize(std::max(it->from, it->to) + 1);

    if(it->fromlabel != -1)
      g[it->from].label = it->fromlabel;
    if(it->tolabel != -1)
      g[it->to].label = it->tolabel;

    g[it->from].push(it->from, it->to, it->elabel);
    if(g.directed == false)
      g[it->to].push(it->to, it->from, it->elabel);
  }

  g.buildEdge();

  return (true);
}

unsigned int
DFSCode::nodeCount(void)
{
  unsigned int nodecount = 0;

  for(DFSCode::iterator it = begin(); it != end(); ++it)
    nodecount = std::max(nodecount, (unsigned int)(std::max(it->from, it->to) + 1));

  return (nodecount);
}



size_t DFSCode::get_serialized_size() const
{
  if(empty()) return sizeof(int);

  return size() * at(0).get_serialized_size() + sizeof(int);
}


size_t DFSCode::get_serialized_size(char *buffer, size_t buffer_size) const
{
  if(buffer_size < sizeof(int)) throw runtime_error("Buffer too small.");
  return *((int *) buffer);
}


size_t DFSCode::serialize(char *buffer, size_t buffer_size) const
{
  if(buffer_size < get_serialized_size()) throw runtime_error("Buffer too small.");

  char *buf = buffer;
  size_t stored = 0;

  // store dfs code element count
  *((int*)(buf + stored)) = size();
  stored += sizeof(int);

  // store each dfs element
  for(int i = 0; i < size(); i++) {
    size_t tmp = at(i).serialize(buf + stored, buffer_size - stored);
    stored += tmp;
  } // for i

  return stored;
} // DFSCode::serialize


size_t DFSCode::deserialize(char *buffer, size_t buffer_size)
{
  if(get_serialized_size(buffer, buffer_size) == 0) return sizeof(int);
  clear();

  int elements = *((int*)buffer);
  size_t read = sizeof(int);

  for(int i = 0; i < elements; i++) {
    DFS d;
    size_t tmp = d.deserialize(buffer + read, buffer_size - read);
    read += tmp;
    push_back(d);
  } // for i
  return read;
} // DFSCode::deserialize


size_t DFSCode::serialize(std::ostream &) const
{
  throw method_unimplemented("DFSCode::serialize");
} // DFSCode::serialize


size_t DFSCode::deserialize(std::istream &)
{
  throw method_unimplemented("DFSCode::deserialize");
} // DFSCode::deserialize



void DFSCode::remove_negative_ones()
{
  if(size() < 1) return;

  int last_vid = at(0).to;

  for(int i = 1; i < size(); i++) {
    if(at(i).from == -1) at(i).from = last_vid;
    if(at(i).is_forward()) last_vid = at(i).to;
  }
} // DFSCode::remove_negative_ones

} // namespace types


