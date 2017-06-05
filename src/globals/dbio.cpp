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

#include <dbio.hpp>
#include <sys/stat.h>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <numeric>
#include <iostream>
#include <logger.hpp>


using std::string;
using std::stringstream;
using std::runtime_error;
using std::find;
using std::make_pair;

using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;
using std::flush;
using std::streampos;

using namespace std;

namespace dbio {


FILE_TYPE ftype2str(const std::string &ftype)
{
  if(ftype == "-txt" || ftype == "txt") return TXT;
  if(ftype == "-adj" || ftype == "adj") return ADJ;
  if(ftype == "-adjp" || ftype == "adjp") return ADJP;
  if(ftype == "-ladjp" || ftype == "ladjp") return LADJP;

  throw runtime_error("unknown filetype");
}

void read_graph(const FILE_TYPE &ftype, const std::string &filename, types::Graph &g)
{
  switch(ftype) {
  case TXT:
    read_graph_txt(filename, g);
    break;
  case ADJ:
    read_graph_adj(filename, g);
    break;
  case LADJP:
    read_graph_local_adjp(filename, g);
    break;
  default:
    throw runtime_error("invalid file type for read_graph");
  } // switch
} // read_graph


void read_graph_txt(const string &filename, types::Graph &g)
{
  ifstream in;
  in.open(filename.c_str(), ios::in);
  try {
    g.read_txt(in);
  } catch(std::runtime_error &e) {
    in.close();
    throw;
  }
  in.close();

} // read_graph_txt


void read_graph_adj(const string &filename, types::Graph &g)
{
  ifstream in;
  in.open(filename.c_str(), ios::in);
  try {
    g.read_adj(in);
  } catch(std::runtime_error &e) {
    in.close();
    throw;
  }
  in.close();
} // read_graph_adj

void read_graph_local_adjp(const string &filename, types::Graph &g)
{
  ifstream in;
  in.open(filename.c_str(), ios::in);
  try {
    g.has_ext_neighbor = false;
    g.read_local_adjp(in);
  } catch(std::runtime_error &e) {
    in.close();
    throw;
  }
  in.close();
} // read_graph_local_adjp


void read_graph_adj_par(const string &graph_filename, const string &partition_filename, int rank, types::Graph &g)
{
  ifstream in_g, in_p;

  //first read the partition ids for all vertices
  //and populate the relevant ones
  in_p.open(partition_filename.c_str(), ios::in);
  try {
    g.has_ext_neighbor = false;
    g.read_partition_info(in_p, rank);
  } catch(std::runtime_error &e) {
    in_p.close();
    throw;
  }

  //read the graph (partition)
  in_g.open(graph_filename.c_str(), ios::in);
  try {
    g.read_adj_par(in_g);
  } catch(std::runtime_error &e) {
    in_g.close();
    throw;
  }
  in_g.close();

  // get the partition ids for the boundary vertices
  in_p.clear();
  in_p.seekg(0,std::ios_base::beg);
  try {
    g.read_partition_info_non_locals(in_p);
  } catch(std::runtime_error &e) {
    in_p.close();
    throw;
  }
  in_p.close();
} // read_graph_adj_par

void write_local_adjp(std::string dirname_out, int rank, types::Graph &g){

  char buf[512];
  mkdir(dirname_out.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);

  std::sprintf(buf, "%s/%d", dirname_out.c_str(),rank);

  std::ofstream out;
  out.open(buf);

  try {
    g.write_local_adjp(out);
  } catch(std::runtime_error &e) {
    out.close();
    throw;
  }

  out.close();
}

void remove_vmap_file(int rank, int thread_id){
  char buf[512];
  std::sprintf(buf, "vmaps_p%dt%d",rank, thread_id);
  std::remove(buf);
}

void remove_vmap_file(const char *filename){
  std::remove(filename);
}

void write_vertex_mappings(int rank, int thread_id, std::vector<std::set<int> > &vmaps, int &position, int &length){

  char* fname = new char[512];
  std::sprintf(fname, "vmaps_p%dt%d",rank, thread_id);
  return write_vertex_mappings(fname, vmaps, position, length);
  delete [] fname;
}

void write_vertex_mappings(const char* filename, std::vector<std::set<int> > &vmaps, int &position, int &length){

  //read mappings in memory buf
  int len = 1 + vmaps.size();
  for(int i = 0; i< vmaps.size(); i++) {
    len += vmaps[i].size();
  }

  int* buf = new int[len];
  int offset = 0;
  buf[offset++] = vmaps.size();

  for(int i = 0; i< vmaps.size(); i++) {
    buf[offset++] = vmaps[i].size();
  }

  for(int i = 0; i< vmaps.size(); i++) {
    for(std::set<int>::iterator it = vmaps[i].begin(); it != vmaps[i].end(); it++) {
      buf[offset++] = *it;
    }
  }

  //write mappings into binary file
  ofstream file (filename, ios::out | ios::binary | ios::app);
  if (file.is_open())
  {
    position = file.tellp();

    file.write ((char*) buf, len * sizeof(int));
    file.close();
  }

  length = len;
  delete [] buf;
}

void read_vertex_mappings(int rank, int thread_id, int position, int length, std::vector<std::set<int> > &vmaps){
  char* fname = new char[512];
  std::sprintf(fname, "vmaps_p%dt%d",rank, thread_id);
  read_vertex_mappings(fname, position, length, vmaps);
  delete [] fname;
}


void read_vertex_mappings(const char* filename, int pos, int len, std::vector<std::set<int> > &vmaps){

  //read mappings in memory buf
  int* buf = new int[len];
  ifstream file(filename, ios::in | ios::binary);
  if (file.is_open())
  {
    file.seekg(pos, ios::beg);
    file.read((char*) buf, len * sizeof(int));
    file.close();
  }

  int offset = 0;
  int vmap_size = buf[offset++];

  vmaps.clear();
  std::vector<int> lengths;
  for(int i = 0; i< vmap_size; i++) {
    std::set<int> t;
    vmaps.push_back(t);
    lengths.push_back(buf[offset++]);
  }

  for(int i = 0; i<vmap_size; i++) {
    for(int j = 0; j<lengths[i]; j++) {
      vmaps[i].insert(buf[offset++]);
    }
  }

  delete [] buf;
}

} // namespace dbio
