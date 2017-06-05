/*
 * This source is part of the single graph mining algorithm.
 *
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



#ifndef __DBIO_HPP__
#define __DBIO_HPP__

#include <fstream>
#include <graph_types.hpp>

namespace dbio {
typedef enum { TXT, ADJ, ADJP, LADJP } FILE_TYPE;
FILE_TYPE ftype2str(const std::string &ftype);

// reads single graph
void read_graph(const FILE_TYPE &ftype, const std::string &filename, types::Graph &g);
void read_graph_txt(const std::string &filename, types::Graph &g);
void read_graph_adj(const std::string &filename, types::Graph &g);
void read_graph_local_adjp(const std::string &filename, types::Graph &g);
void read_graph_adj_par(const std::string &graph_filename, const std::string &partition_filename, int rank, types::Graph &g);

//void write_graph_adj_par(const std::string &output_filename, int rank, types::Graph &g);
void write_local_adjp(std::string dirname_out, int rank, types::Graph &g);

//vmap files
void remove_vmap_file(int rank, int thread_id);
void remove_vmap_file(const char *filename);
void write_vertex_mappings(int rank, int thread_id, std::vector<std::set<int> > &vmaps, int &position, int &length);
void write_vertex_mappings(const char* filename, std::vector<std::set<int> > &vmaps, int &position, int &length);
void read_vertex_mappings(int rank, int thread_id, int pos, int len, std::vector<std::set<int> > &vmaps);
void read_vertex_mappings(const char* filename, int pos, int len, std::vector<std::set<int> > &vmaps);

} // namespace dbio

#endif
