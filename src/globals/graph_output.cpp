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

#include <graph_output.hpp>
#include <sstream>
#include <fstream>

using std::stringstream;

graph_output::~graph_output()
{
}

void graph_counting_output::output_graph(const types::Graph &g, size_t support) {
  counter += 1.0L;
} // output_sequence


void graph_counting_output::output_graph(const types::DFSCode &code, size_t support)
{
  counter += 1.0L;
}

void graph_counting_output::output_graph(const types::DFSCode &code)
{
  counter += 1.0L;
}

graph_counting_output::graph_counting_output() {
  counter = 0;
}

std::string graph_counting_output::to_string() {
  std::stringstream ss;
  ss << counter;
  return ss.str();
} // to_string


double graph_counting_output::get_count()
{
  return counter;
}


graph_file_output::graph_file_output(const std::string &fname) {
  filename = fname;
  out = new std::ofstream();
  ((std::ofstream*)out)->open(filename.c_str(), std::ios::out | std::ios::trunc);
}


void graph_file_output::output_graph(const types::Graph &g, size_t support)
{
  types::DFSCode dfs = g.get_min_dfs_code();
  output_graph(dfs);
  //(*out) << g.to_string().c_str() << std::endl;
  //count += 1.0L;
}

void graph_file_output::output_graph(const types::DFSCode &code, size_t support)
{
  (*out) << code.to_string(false) << ":" << support << std::endl;
  count += 1.0L;
}

void graph_file_output::output_graph(const types::DFSCode &code)
{
  (*out) << code.to_string(false) << std::endl;
  count += 1.0L;
}

double graph_file_output::get_count()
{
  return count;
}

std::string graph_file_output::to_string()
{
  stringstream ss;
  ss << "graph_file_output, filename: " << filename;
  return ss.str();
}


