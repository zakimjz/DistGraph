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

#ifndef __DFS_CODE_HPP__
#define __DFS_CODE_HPP__

#include <types.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <serializable.hpp>

namespace types {

class Graph;
class Projected;
class History;

class DFS {
public:
  int from;
  int to;
  int fromlabel;
  int elabel;
  int tolabel;


  friend bool operator==(const DFS &d1, const DFS &d2) {
    return (d1.from == d2.from && d1.to == d2.to && d1.fromlabel == d2.fromlabel
            && d1.elabel == d2.elabel && d1.tolabel == d2.tolabel);
  }
  friend bool operator!=(const DFS &d1, const DFS &d2) {
    return (!(d1 == d2));
  }

  friend std::ostream &operator<<(std::ostream &out, const DFS &d) {
    out << d.to_string().c_str();
    return out;
  }

  friend bool operator<(const DFS &d1, const DFS &d2){
    if(d1.from < d2.from) return true;
    if(d1.from > d2.from) return false;

    if(d1.to < d2.to) return true;
    if(d1.to > d2.to) return false;

    if(d1.fromlabel < d2.fromlabel) return true;
    if(d1.fromlabel > d2.fromlabel) return false;

    if(d1.elabel < d2.elabel) return true;
    if(d1.elabel > d2.elabel) return false;

    if(d1.tolabel < d2.tolabel) return true;
    if(d1.tolabel > d2.tolabel) return false;

    return false;
  }

  DFS() : from(0), to(0), fromlabel(0), elabel(0), tolabel(0) {
  }
  DFS(int from, int to, int fromlabel, int elabel, int tolabel) : from(from), to(to), fromlabel(fromlabel), elabel(elabel), tolabel(tolabel) {
  }
  DFS(char *buffer, int size);

  DFS(const DFS &d) : from(d.from), to(d.to), fromlabel(d.fromlabel), elabel(d.elabel), tolabel(d.tolabel) {
  }

  std::string to_string(bool print_edge_type = true) const;

  bool is_forward() const;
  bool is_backward() const;

  static size_t deserialize(DFS &result, char *buf, size_t bufsize);
  static size_t serialize(DFS &result, char *buf, size_t bufsize);
  static size_t get_serialized_dfs_code_size(DFS &result);
  static DFS parse_from_string(const char *str_dfs_code);


  virtual size_t get_serialized_size() const;
  virtual size_t get_serialized_size(char *buffer, size_t buffer_size) const;
  virtual size_t serialize(char *buffer, size_t buffer_size) const;
  virtual size_t deserialize(char *buffer, size_t buffer_size);


  virtual size_t serialize(std::ostream &) const;
  virtual size_t deserialize(std::istream &);
};

struct DFS_less_then {
  bool bckwrd_bckwrd_less(const DFS &d1, const DFS &d2) const;
  bool frwrd_bckwrd_less(const DFS &d1, const DFS &d2) const;
  bool frwrd_frwrd_less(const DFS &d1, const DFS &d2) const;

  bool operator()(const DFS &d1, const DFS &d2) const;
};


struct DFS_less_then_fast {
  bool operator()(const DFS &d1, const DFS &d2) const {
    if(d1.from < d2.from) return true;
    if(d1.from > d2.from) return false;

    if(d1.to < d2.to) return true;
    if(d1.to > d2.to) return false;

    if(d1.fromlabel < d2.fromlabel) return true;
    if(d1.fromlabel > d2.fromlabel) return false;

    if(d1.elabel < d2.elabel) return true;
    if(d1.elabel > d2.elabel) return false;

    if(d1.tolabel < d2.tolabel) return true;
    if(d1.tolabel > d2.tolabel) return false;


    return false;
  }
};


/**
 * Standard equal, as defined in DFS. The only difference is that it
 * is functor.
 */
struct DFS_std_equal {
  bool operator()(const DFS &d1, const DFS &d2) const {
    return d1 == d2;
  }
};


struct DFS_std_not_equal {
  bool operator()(const DFS &d1, const DFS &d2) const {
    return d1 != d2;
  }
};

/**
 * this is a special version of the == operator that compares the DFS
 * structure only partially, depending on whether it is forward or
 * backward edge.
 */
struct DFS_partial_equal {
  bool operator()(const DFS &d1, const DFS &d2) const {
    if(d1.from != d2.from || d1.to != d2.to) return false;
    if(d1.fromlabel != -1 && d2.fromlabel != -1 && d1.fromlabel != d2.fromlabel) return false;
    if(d1.tolabel != -1 && d2.tolabel != -1 && d1.tolabel != d2.tolabel) return false;
    if(d1.elabel != d2.elabel) return false;
    return true;
  } // operator()
};


/**
 * this is a special version of the == operator that compares the DFS
 * structure only partially, depending on whether it is forward or
 * backward edge. INTERNALLY USES DFS_equal.
 */
struct DFS_partial_not_equal {
  DFS_partial_equal eq;
  bool operator()(const DFS &d1, const DFS &d2) const {
    return !eq(d1, d2);
  }
};

class DFSCode;
std::ostream &operator<<(std::ostream &out, const DFSCode &code);

struct DFSCode : public std::vector<DFS>, public serializable {
private:
  // right-most path
  RMPath rmpath;
public:
  const RMPath &get_rmpath() const {
    return rmpath;
  }
  const RMPath &buildRMPath();

  // Convert current DFS code into a graph.
  bool toGraph(Graph &) const;

  // Clear current DFS code and build code from the given graph.
  //void fromGraph(Graph &g);

  // Return number of nodes in the graph.
  unsigned int nodeCount(void);

  DFSCode &operator=(const DFSCode &other) {
    if(this == &other) return *this;
    std::vector<DFS>::operator=(other);
    rmpath = other.rmpath;
    return *this;
  }

  friend bool operator==(const DFSCode &d1, const DFSCode &d2) {
    if(d1.size() != d2.size())
      return false;
    for(int i = 0; i < d1.size(); i++) {
      if(d1[i] != d2[i])
        return false;
    }
    return true;
  }

  friend bool operator<(const DFSCode &d1, const DFSCode &d2) {
    if(d1.size() < d2.size())
      return true;
    else if(d1.size() > d2.size())
      return false;

    for(int i = 0; i < d1.size(); i++) {
      if(d1[i] < d2[i])
        return true;
      else if(d2[i] < d1[i])
        return false;
    }
    return false;         //equal
  }

  friend std::ostream &operator<<(std::ostream &out, const DFSCode &code);

  void push(int from, int to, int fromlabel, int elabel, int tolabel) {
    resize(size() + 1);
    DFS &d = (*this)[size() - 1];

    d.from = from;
    d.to = to;
    d.fromlabel = fromlabel;
    d.elabel = elabel;
    d.tolabel = tolabel;
  }
  void pop() {
    resize(size() - 1);
  }
  std::ostream &write(std::ostream &) const;  // write
  std::string to_string(bool print_edge_type = true) const;

  bool dfs_code_is_min() const;

  static types::DFSCode read_from_str(const std::string &str);

  virtual size_t get_serialized_size() const;
  virtual size_t get_serialized_size(char *buffer, size_t buffer_size) const;
  virtual size_t serialize(char *buffer, size_t buffer_size) const;
  virtual size_t deserialize(char *buffer, size_t buffer_size);

  virtual size_t serialize(std::ostream &) const;
  virtual size_t deserialize(std::istream &);

  void remove_negative_ones();

protected:
  static bool dfs_code_is_min_internal(types::Projected &projected, const DFSCode &dfs_code, DFSCode &min_dfs_code, types::Graph &graph_dfs_code);
};

} // namespace types

#endif

