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

#ifndef __GRAPH_OUTPUT_HPP__
#define __GRAPH_OUTPUT_HPP__


#include <graph_types.hpp>
#include <fstream>


class graph_output {
public:
  virtual void output_graph(const types::Graph &g, size_t support) = 0;
  virtual void output_graph(const types::DFSCode &code, size_t support) = 0;
  virtual void output_graph(const types::DFSCode &code) = 0;
  virtual double get_count() = 0;
  virtual std::string to_string() = 0;
  virtual ~graph_output();
};


class graph_counting_output : public graph_output {
  double counter;
public:
  graph_counting_output();
  virtual void output_graph(const types::Graph &g, size_t support);
  virtual void output_graph(const types::DFSCode &code, size_t support);
  virtual void output_graph(const types::DFSCode &code);
  virtual double get_count();
  virtual std::string to_string();
};


class graph_file_output : public graph_output {
  std::ostream *out;
  std::string filename;
  double count;
public:
  graph_file_output(const std::string &fname);
  virtual void output_graph(const types::Graph &g, size_t support);
  virtual void output_graph(const types::DFSCode &code, size_t support);
  virtual void output_graph(const types::DFSCode &code);
  virtual double get_count();
  virtual std::string to_string();
};

#endif
