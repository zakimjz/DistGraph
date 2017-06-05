/*
 * This source is part of the single graph mining algorithm.
 *
 * Copyright 2014-2016 Nilothpal Talukder
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


#include <global_support_handler.hpp>
#include <unistd.h>
#include <string>
#include <cmath>
#include <algorithm>
#include <alloc_tools.hpp>
#include <dbio.hpp>

using namespace std;

namespace algs {

void global_support_handler::init_global_support_handler(unsigned int rank, unsigned int numtasks, unsigned int numthreads)
{
  this->rank = rank;
  this->numtasks = numtasks;
  this->numthreads = numthreads;
  this->out = 0;
  this->num_partitions = numtasks;

  hash_sendrecv_time = 0;
  compute_support_time = 0;
  unfinished_sendrecv_time = 0;
  global_pattern_sendrecv_time = 0;

  for(int i = 0; i< numthreads; i++) {
    std::map<string, std::vector<int> > t;
    external_supports_threads.push_back(t);
    internal_supports_threads.push_back(t);
    supports_locally_frequent_threads.push_back(t);
    supports_locally_not_frequent_threads.push_back(t);
  }

  frequent_patterns = 0;
  char prefix[1000];
  sprintf(prefix, "vmaps-p%d-t0",rank);
  vmap_filename = std::string(prefix);

  logger = Logger::get_logger("GLOBAL_SUPPORT_HANDLER");
  TRACE(*logger, "global support handler initialized");
}

void global_support_handler::set_min_support(int minsupp){
  this->minsupp = minsupp;
}

void global_support_handler::set_num_partitions(int num_partitions){
  this->num_partitions = num_partitions;
}

void global_support_handler::set_graph_output(graph_output * gout)
{

  DEBUG(*logger, " inside set graph output");
  out = gout;
}

void global_support_handler::set_vmap_filename(std::string filename){
  vmap_filename = filename;
}

void global_support_handler::set_optimization_flags(bool global_support_optimization){
  this->global_support_optimization = global_support_optimization;
}

void global_support_handler::insert_local_support_data(types::DFSCode DFS_CODE, std::vector<int> vertex_supp, bool is_locally_frequent){
  string dfscode = DFS_CODE.to_string();

  DEBUG(*logger, " inserting local support data for " << DFS_CODE.to_string() << " size = " << vertex_supp.size());
  omp_set_lock(&lock_gs);
  if(is_locally_frequent) {
    supports_locally_frequent[dfscode] = vertex_supp;
  } else {
    supports_locally_not_frequent[dfscode] = vertex_supp;
  }
  omp_unset_lock(&lock_gs);
}

void global_support_handler::insert_local_support_data_no_lock(types::DFSCode DFS_CODE, std::vector<int> vertex_supp, bool is_locally_frequent){
  string dfscode = DFS_CODE.to_string();
  int thread_id = omp_get_thread_num();

  DEBUG(*logger, " inserting local support data for " << DFS_CODE.to_string() << " size = " << vertex_supp.size());
  if(is_locally_frequent) {
    supports_locally_frequent_threads[thread_id][dfscode] = vertex_supp;
  } else {
    supports_locally_not_frequent_threads[thread_id][dfscode] = vertex_supp;
  }
}

int global_support_handler::compute_local_support_and_store_mappings(types::Embeddings &embeddings, types::DFSCode DFS_CODE, bool store_in_file, std::vector<int> &vertex_support){
  string dfscode = DFS_CODE.to_string();

  int maxtoc = 0;
  for(int i = 0; i < DFS_CODE.size(); i++) {
    if(DFS_CODE[i].is_forward())
      maxtoc++;
  }

  std::vector<std::set<int> > maps(maxtoc + 1);

  //iterated through the all the embeddings
  for(types::Embeddings::iterator cur = embeddings.begin(); cur != embeddings.end(); ++cur) {

    types::Emb *em = &(*cur);
    int dfsindex = DFS_CODE.size() - 1;
    while(em != NULL) {
      if(DFS_CODE[dfsindex].to > DFS_CODE[dfsindex].from) {       //forward edge
        maps[DFS_CODE[dfsindex].to].insert((*em).edge.to);
      }
      if(!em->prev) {
        maps[DFS_CODE[dfsindex].from].insert((*em).edge.from);
      }

      //if(maxtoc < DFS_CODE[dfsindex].to)
      //maxtoc = DFS_CODE[dfsindex].to;

      em = em->prev;
      dfsindex--;
    }
  }

  int max = 0;
  for(int i = 0; i<=maxtoc; i++) {
    vertex_support.push_back(maps[i].size());
    if(maps[i].size() > max)
      max = maps[i].size();
  }


  int pos, len;
  omp_set_lock(&lock_gs);
  if(store_in_file) {
    //dbio::write_vertex_mappings(rank,0,maps, pos, len);
    dbio::write_vertex_mappings(vmap_filename.c_str(),maps, pos, len);
    vertex_mappings_file_pos[dfscode] = pos;
    vertex_mappings_size_in_file[dfscode] = len;
  } else {
    vertex_mappings[dfscode] = maps;
  }

  omp_unset_lock(&lock_gs);

  return max;
}

int global_support_handler::compute_local_support(types::Embeddings &embeddings, types::DFSCode DFS_CODE, std::vector<int> &vertex_support){
  string dfscode = DFS_CODE.to_string();

  int maxtoc = 0;
  for(int i = 0; i < DFS_CODE.size(); i++) {
    if(DFS_CODE[i].is_forward())
      maxtoc++;
  }

  std::vector<std::set<int> > maps(maxtoc + 1);

  //iterated through the all the embeddings
  for(types::Embeddings::iterator cur = embeddings.begin(); cur != embeddings.end(); ++cur) {

    types::Emb *em = &(*cur);
    int dfsindex = DFS_CODE.size() - 1;
    while(em != NULL) {
      if(DFS_CODE[dfsindex].to > DFS_CODE[dfsindex].from) {       //forward edge
        maps[DFS_CODE[dfsindex].to].insert((*em).edge.to);
      }
      if(!em->prev) {
        maps[DFS_CODE[dfsindex].from].insert((*em).edge.from);
      }

      em = em->prev;
      dfsindex--;
    }
  }

  int max = 0;
  for(int i = 0; i<=maxtoc; i++) {
    vertex_support.push_back(maps[i].size());
    if(maps[i].size() > max)
      max = maps[i].size();
  }


  int pos, len;
  omp_set_lock(&lock_gs);
  vertex_mappings[dfscode] = maps;
  omp_unset_lock(&lock_gs);

  return max;
}


int global_support_handler::compute_local_support_and_store_external_mappings(types::Graph &graph, types::Embeddings &embeddings, types::DFSCode DFS_CODE, std::vector<int> &vertex_support){
  string dfscode = DFS_CODE.to_string();

  int maxtoc = 0;
  for(int i = 0; i < DFS_CODE.size(); i++) {
    if(DFS_CODE[i].is_forward())
      maxtoc++;
  }

  std::vector<std::set<int> > maps(maxtoc + 1);
  std::vector<int> int_supports(maxtoc + 1), ext_supports(maxtoc + 1);
  //std::vector<std::map<int, std::set<int> > > external_maps_per_part(maxtoc + 1);
  std::vector<std::set<int> > external_maps(maxtoc + 1);

  //iterated through the all the embeddings
  for(types::Embeddings::iterator cur = embeddings.begin(); cur != embeddings.end(); ++cur) {

    types::Emb *em = &(*cur);
    int dfsindex = DFS_CODE.size() - 1;
    while(em != NULL) {
      if(DFS_CODE[dfsindex].to > DFS_CODE[dfsindex].from) {       //forward edge
        maps[DFS_CODE[dfsindex].to].insert((*em).edge.to);
        if(graph[(*em).edge.to].orig_part_id != rank % num_partitions) {
          external_maps[DFS_CODE[dfsindex].to].insert((*em).edge.to);
        }
      }
      if(!em->prev) {
        maps[DFS_CODE[dfsindex].from].insert((*em).edge.from);
        if(graph[(*em).edge.from].orig_part_id != rank % num_partitions) {
          external_maps[DFS_CODE[dfsindex].from].insert((*em).edge.from);
        }
      }

      em = em->prev;
      dfsindex--;
    }
  }

  int max = 0;
  for(int i = 0; i<=maxtoc; i++) {
    vertex_support.push_back(maps[i].size());
    ext_supports[i] = external_maps[i].size();
    int_supports[i] = maps[i].size() - ext_supports[i];
    if(maps[i].size() > max)
      max = maps[i].size();
  }

  int pos, len;
  //std::vector<std::set<int> > maps_r;
  omp_set_lock(&lock_gs);
  //vertex_mappings[dfscode] = maps;
  external_mappings[dfscode] = external_maps;
  //external_mappings_per_part[dfscode] = external_maps_per_part;
  internal_supports[dfscode] = int_supports;
  external_supports[dfscode] = ext_supports;
  //dbio::write_vertex_mappings(rank,0,maps, pos, len);
  //vertex_mappings_file_pos[dfscode] = pos;
  //vertex_mappings_size_in_file[dfscode] = len;
  //dbio::read_vertex_mappings(rank, 0, pos, len, maps_r);
  DEBUG(*logger, " internal supports " << utils::print_vector(internal_supports[dfscode]) << " external supports " << utils::print_vector(external_supports[dfscode]));
  omp_unset_lock(&lock_gs);


  /*if (maps.size() != maps_r.size())
          CRITICAL_ERROR(*logger,"SIZES don't match"<< maps.size()<<" "<<maps_r.size());
     for(int i = 0; i < maps.size(); i++ )
          for(std::set<int>::iterator it = maps[i].begin(), it2 = maps_r[i].begin(); it != maps[i].end(); it++, it2++)
                  if(*it != *it2)
                          CRITICAL_ERROR(*logger,"Vertex ids don't match"<< *it<<" "<<*it2);
   */

  return max;
}

int global_support_handler::compute_local_support_and_store_int_ext_supports(types::Graph &graph, types::Embeddings &embeddings, types::DFSCode DFS_CODE, std::vector<int> &vertex_support){
  string dfscode = DFS_CODE.to_string();

  int maxtoc = 0;
  for(int i = 0; i < DFS_CODE.size(); i++) {
    if(DFS_CODE[i].is_forward())
      maxtoc++;
  }

  std::vector<std::set<int> > maps(maxtoc + 1);
  std::vector<std::set<int> > external_maps(maxtoc + 1);
  std::vector<int> int_supports, ext_supports;

  //iterated through the all the embeddings
  for(types::Embeddings::iterator cur = embeddings.begin(); cur != embeddings.end(); ++cur) {

    types::Emb *em = &(*cur);
    int dfsindex = DFS_CODE.size() - 1;
    while(em != NULL) {
      if(DFS_CODE[dfsindex].to > DFS_CODE[dfsindex].from) {       //forward edge
        maps[DFS_CODE[dfsindex].to].insert((*em).edge.to);
        if(graph[(*em).edge.to].orig_part_id != rank % num_partitions)
          external_maps[DFS_CODE[dfsindex].to].insert((*em).edge.to);
      }
      if(!em->prev) {
        maps[DFS_CODE[dfsindex].from].insert((*em).edge.from);
        if(graph[(*em).edge.from].orig_part_id != rank % num_partitions)
          external_maps[DFS_CODE[dfsindex].from].insert((*em).edge.from);
      }

      em = em->prev;
      dfsindex--;
    }
  }

  int max = 0;
  for(int i = 0; i<=maxtoc; i++) {
    vertex_support.push_back(maps[i].size());
    ext_supports.push_back(external_maps[i].size());
    int_supports.push_back(maps[i].size() - external_maps[i].size());
    if(maps[i].size() > max)
      max = maps[i].size();
  }

  int pos, len;
  omp_set_lock(&lock_gs);
  internal_supports[dfscode] = int_supports;
  external_supports[dfscode] = ext_supports;
  //dbio::write_vertex_mappings(rank,0,maps, pos, len);
  //vertex_mappings_file_pos[dfscode] = pos;
  //vertex_mappings_size_in_file[dfscode] = len;
  omp_unset_lock(&lock_gs);
  return max;
}

int global_support_handler::compute_local_support_and_store_int_ext_supports_no_lock(types::Graph &graph, types::Embeddings &embeddings, types::DFSCode DFS_CODE, std::vector<int> &vertex_support){
  string dfscode = DFS_CODE.to_string();
  int thread_id = omp_get_thread_num();

  int maxtoc = 0;
  for(int i = 0; i < DFS_CODE.size(); i++) {
    if(DFS_CODE[i].is_forward())
      maxtoc++;
  }

  std::vector<std::set<int> > maps(maxtoc + 1);
  std::vector<std::set<int> > external_maps(maxtoc + 1);
  std::vector<int> int_supports, ext_supports;

  //iterated through the all the embeddings
  for(types::Embeddings::iterator cur = embeddings.begin(); cur != embeddings.end(); ++cur) {

    types::Emb *em = &(*cur);
    int dfsindex = DFS_CODE.size() - 1;
    while(em != NULL) {
      if(DFS_CODE[dfsindex].to > DFS_CODE[dfsindex].from) {       //forward edge
        maps[DFS_CODE[dfsindex].to].insert((*em).edge.to);
        if(graph[(*em).edge.to].orig_part_id != rank % num_partitions)
          external_maps[DFS_CODE[dfsindex].to].insert((*em).edge.to);
      }
      if(!em->prev) {
        maps[DFS_CODE[dfsindex].from].insert((*em).edge.from);
        if(graph[(*em).edge.from].orig_part_id != rank % num_partitions)
          external_maps[DFS_CODE[dfsindex].from].insert((*em).edge.from);
      }

      em = em->prev;
      dfsindex--;
    }
  }

  int max = 0;
  for(int i = 0; i<=maxtoc; i++) {
    vertex_support.push_back(maps[i].size());
    ext_supports.push_back(external_maps[i].size());
    int_supports.push_back(maps[i].size() - external_maps[i].size());
    if(maps[i].size() > max)
      max = maps[i].size();
  }

  int pos, len;
//omp_set_lock(&lock_gs);
  internal_supports_threads[thread_id][dfscode] = int_supports;
  external_supports_threads[thread_id][dfscode] = ext_supports;
  //dbio::write_vertex_mappings(rank,0,maps, pos, len);
  //vertex_mappings_file_pos[dfscode] = pos;
  //vertex_mappings_size_in_file[dfscode] = len;
  //omp_unset_lock(&lock_gs);
  return max;
}


int global_support_handler::compute_local_support_and_store_mappings_int_ext_supports(types::Graph &graph, types::Embeddings &embeddings, types::DFSCode DFS_CODE, bool store_in_file, std::vector<int> &vertex_support){
  string dfscode = DFS_CODE.to_string();

  int maxtoc = 0;
  for(int i = 0; i < DFS_CODE.size(); i++) {
    if(DFS_CODE[i].is_forward())
      maxtoc++;
  }

  std::vector<std::set<int> > maps(maxtoc + 1);
  std::vector<std::set<int> > external_maps(maxtoc + 1);
  std::vector<int> int_supports, ext_supports;

  //iterated through the all the embeddings
  for(types::Embeddings::iterator cur = embeddings.begin(); cur != embeddings.end(); ++cur) {

    types::Emb *em = &(*cur);
    int dfsindex = DFS_CODE.size() - 1;
    while(em != NULL) {
      if(DFS_CODE[dfsindex].to > DFS_CODE[dfsindex].from) {       //forward edge
        maps[DFS_CODE[dfsindex].to].insert((*em).edge.to);
        if(graph[(*em).edge.to].orig_part_id != rank % num_partitions)
          external_maps[DFS_CODE[dfsindex].to].insert((*em).edge.to);
      }
      if(!em->prev) {
        maps[DFS_CODE[dfsindex].from].insert((*em).edge.from);
        if(graph[(*em).edge.from].orig_part_id != rank % num_partitions)
          external_maps[DFS_CODE[dfsindex].from].insert((*em).edge.from);
      }

      em = em->prev;
      dfsindex--;
    }
  }

  int max = 0;
  for(int i = 0; i<=maxtoc; i++) {
    vertex_support.push_back(maps[i].size());
    ext_supports.push_back(external_maps[i].size());
    int_supports.push_back(maps[i].size() - external_maps[i].size());
    if(maps[i].size() > max)
      max = maps[i].size();
  }

  int pos, len;
  omp_set_lock(&lock_gs);
  internal_supports[dfscode] = int_supports;
  external_supports[dfscode] = ext_supports;
  if(store_in_file) {
    //dbio::write_vertex_mappings(rank,0,maps, pos, len);
    dbio::write_vertex_mappings(vmap_filename.c_str(),maps, pos, len);
    vertex_mappings_file_pos[dfscode] = pos;
    vertex_mappings_size_in_file[dfscode] = len;
  }else{
    vertex_mappings[dfscode] = maps;
  }
  omp_unset_lock(&lock_gs);

  return max;
}

void global_support_handler::combine_threads_map_vectors(){

  for(int i = 0; i< numthreads; i++) {
    supports_locally_frequent.insert(supports_locally_frequent_threads[i].begin(), supports_locally_frequent_threads[i].end());
    supports_locally_not_frequent.insert(supports_locally_not_frequent_threads[i].begin(), supports_locally_not_frequent_threads[i].end());
    external_supports.insert(external_supports_threads[i].begin(), external_supports_threads[i].end());
    internal_supports.insert(internal_supports_threads[i].begin(), internal_supports_threads[i].end());

    supports_locally_frequent_threads[i].clear();
    supports_locally_not_frequent_threads[i].clear();
    external_supports_threads[i].clear();
    internal_supports_threads[i].clear();
  }

}

bool global_support_handler::get_stored_vertex_mappings(types::DFSCode DFS_CODE, bool in_file, std::vector<std::set<int> > &vmap){
  string dfscode = DFS_CODE.to_string();
  bool found = false;
  if(in_file) {
    if(vertex_mappings_file_pos.count(dfscode) > 0) {
      //dbio::read_vertex_mappings(rank, 0, vertex_mappings_file_pos[dfscode], vertex_mappings_size_in_file[dfscode], vmap);
      dbio::read_vertex_mappings(vmap_filename.c_str(), vertex_mappings_file_pos[dfscode], vertex_mappings_size_in_file[dfscode], vmap);
      found = true;
    }
  }else{
    if(vertex_mappings.count(dfscode) > 0) {
      vmap = vertex_mappings[dfscode];
      found = true;
    }
  }
  return found;
}


void global_support_handler::insert_local_support_for_globally_frequent(types::DFSCode DFS_CODE, int sup){
  string dfscode = DFS_CODE.to_string();
  DEBUG(*logger, " inserting local support data for globally frequent " << DFS_CODE.to_string() << " sup = " << sup);
  local_supports_globally_frequent[dfscode] = sup;
}

int global_support_handler::get_local_support_for_globally_frequent(types::DFSCode DFS_CODE){
  string dfscode = DFS_CODE.to_string();
  if(local_supports_globally_frequent.count(dfscode) > 0 )
    return local_supports_globally_frequent[dfscode];
  else
    return -1;
}

void global_support_handler::remove_local_support_data(types::DFSCode DFS_CODE, bool is_locally_frequent){

  string dfscode = DFS_CODE.to_string();
  std::map<string, std::vector<int> >::iterator it;

  omp_set_lock(&lock_gs);
  if(is_locally_frequent) {
    it = supports_locally_frequent.find(dfscode);
    supports_locally_frequent.erase(it);
  } else {
    it = supports_locally_not_frequent.find(dfscode);
    supports_locally_not_frequent.erase(it);
  }
  omp_unset_lock(&lock_gs);
}


void global_support_handler::gather_and_compute_global_supports(int dfscode_size){

  //INFO(*logger, "rank = "<< rank << " level = "<< dfscode_size <<" locally frequent patterns count = " << supports_locally_frequent.size() << " locally not frequent patterns count = " << supports_locally_not_frequent.size());
  this->dfscode_size = dfscode_size;

  //print_average_vertex_mappings_info();

  //timeval time1, time2;
  //gettimeofday(&time1, 0);
  hash_and_sendrecv_locally_frequent_patterns();

  //gettimeofday(&time2, 0);
  //hash_sendrecv_time += utils::get_time_diff(time1, time2);

  compute_globally_frequent_patterns();

  //gettimeofday(&time1, 0);
  //compute_support_time += utils::get_time_diff(time2, time1);

  sendrecv_unfinished_frequent_patterns();

  //gettimeofday(&time2, 0);
  //unfinished_sendrecv_time += utils::get_time_diff(time1, time2);

  compute_globally_frequent_patterns();

  //gettimeofday(&time1, 0);
  //compute_support_time += utils::get_time_diff(time2, time1);

  //broadcast globally frequent patterns
  send_globally_frequent_patterns();

  //gettimeofday(&time2, 0);
  //global_pattern_sendrecv_time += utils::get_time_diff(time1, time2);

  //MPI_Barrier(MPI_COMM_WORLD);
  INFO(*logger, " TOTAL globally frequent dfscodes = " << globally_frequent_dfscodes.size());

  if(num_partitions < numtasks)
    send_globally_frequent_patterns_to_rank_group();

  //clear all data structures
  hashed_dfscodes.clear();
  global_supports.clear();
  globally_frequent_dfscodes_hashed.clear();

  supports_locally_frequent.clear();
  supports_locally_not_frequent.clear();
}

void global_support_handler::gather_estimate_compute_global_supports(int dfscode_size, types::Graph &graph){

  //INFO(*logger, "rank = "<< rank << " level = "<< dfscode_size <<" locally frequent patterns count = " << supports_locally_frequent.size() << " locally not frequent patterns count = " << supports_locally_not_frequent.size());
  this->dfscode_size = dfscode_size;

  //uncomment this to get all locally frequent and non frequent patterns count
  //std::vector<types::DFSCode> a;
  //get_all_locally_not_frequent_patterns(a);


  //print_average_vertex_mappings_info();

  timeval time1, time2;
  gettimeofday(&time1, 0);

  if(!global_support_optimization)
    hash_and_sendrecv_locally_frequent_patterns();
  else
    hash_and_sendrecv_locally_frequent_patterns_with_int_ext_map_cardinalities();

  gettimeofday(&time2, 0);
  hash_sendrecv_time += utils::get_time_diff(time1, time2);

  gettimeofday(&time1, 0);
  //MPI_Barrier(MPI_COMM_WORLD);
  if(!global_support_optimization)
    compute_globally_frequent_patterns();
  else
    compute_globally_frequent_patterns_with_int_ext_map_cardinalities();

  gettimeofday(&time2, 0);
  compute_support_time += utils::get_time_diff(time1, time2);

  gettimeofday(&time1, 0);
  if(!global_support_optimization)
    sendrecv_unfinished_frequent_patterns();
  else
    sendrecv_unfinished_frequent_patterns_with_int_ext_map_cardinalities();
  //sendrecv_unfinished_frequent_patterns_with_int_ext_map_cardinalities2();

  gettimeofday(&time2, 0);
  unfinished_sendrecv_time += utils::get_time_diff(time1, time2);

  gettimeofday(&time1, 0);
  if(!global_support_optimization)
    compute_globally_frequent_patterns();
  else
    compute_globally_frequent_patterns_with_int_ext_map_cardinalities();

  gettimeofday(&time2, 0);
  compute_support_time += utils::get_time_diff(time1, time2);

  gettimeofday(&time1, 0);
  //broadcast globally frequent patterns
  if(!global_support_optimization)
    send_globally_frequent_patterns();
  else
    send_globally_frequent_and_estimated_patterns();

  gettimeofday(&time2, 0);
  global_pattern_sendrecv_time += utils::get_time_diff(time1, time2);

  INFO(*logger, " TOTAL globally frequent dfscodes = " << globally_frequent_dfscodes.size());

  if(num_partitions < numtasks)
    send_globally_frequent_patterns_to_rank_group();


  //clear all data structures
  hashed_dfscodes.clear();
  global_supports.clear();
  global_supports_internal.clear();
  global_supports_external.clear();
  global_supports_external_max.clear();
  globally_frequent_dfscodes_hashed.clear();
  external_mappings.clear();
  marked_actual_globally_frequent_patterns.clear();

  supports_locally_frequent.clear();
  supports_locally_not_frequent.clear();
  internal_supports.clear();
  external_supports.clear();

}

void global_support_handler::get_stored_dfscode_external_vertex_mappings(std::map<string, std::vector<std::vector<std::set<int> > > > &vertex_mappings){
  vertex_mappings = stored_dfscode_external_vertex_mappings;
}

void global_support_handler::clear_stored_dfscode_external_vertex_mappings(){
  stored_dfscode_external_vertex_mappings.clear();
}



void global_support_handler::gather_and_compute_global_supports_with_mappings(int dfscode_size, types::Graph &graph){

  //INFO(*logger, "rank = "<< rank << " level = "<< dfscode_size <<" locally frequent patterns count = " << supports_locally_frequent.size() << " locally not frequent patterns count = " << supports_locally_not_frequent.size());
  this->dfscode_size = dfscode_size;

  compute_globally_frequent_patterns_with_mappings(graph);

  INFO(*logger, " TOTAL globally frequent dfscodes = " << globally_frequent_dfscodes.size());

  //clear all data structures
  global_supports.clear();
  supports_locally_frequent.clear();
  supports_locally_not_frequent.clear();
  vertex_mappings.clear();
}

void global_support_handler::hash_and_sendrecv_locally_frequent_patterns(){

  // hash locally frequent dfs codes
  for(std::map<string, std::vector<int> >::iterator it = supports_locally_frequent.begin();
      it != supports_locally_frequent.end(); ++it) {

    int hashed_rank = utils::hash(it->first) % numtasks;
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);

    if(hashed_rank == rank) {            //belongs to me? initialize global support
      global_supports[DFS_CODE.to_string()] = it->second;
      /*for(int i = 0; i < global_supports[it->first].size(); i++){
              INFO(*logger, " DFSCODE "<<it->first<< " i = " << i << " total = "<< global_supports[it->first][i]);
         }*/
      continue;
    }

    if(hashed_dfscodes.count(hashed_rank) == 0) {
      std::vector<types::DFSCode> t;
      hashed_dfscodes[hashed_rank] = t;
    }

    hashed_dfscodes[hashed_rank].push_back(DFS_CODE);
  }

  int *sendcounts = new int[numtasks];
  int *recvcounts = new int[numtasks];
  int total_sendcodes = 0;
  int total_recvcodes = 0;

  for(int i = 0; i<numtasks; i++) {
    if(rank == i || hashed_dfscodes.count(i) == 0)
      sendcounts[i] = 0;
    else
      sendcounts[i] = hashed_dfscodes[i].size();
    total_sendcodes += sendcounts[i];
    recvcounts[i] = 0;
  }

  alltoall_msg(sendcounts, recvcounts, 1);
  DEBUG(*logger, "sent counts " << utils::print_array(sendcounts, numtasks) << " recv counts " << utils::print_array(recvcounts, numtasks));

  for(int i = 0; i<numtasks; i++)
    total_recvcodes += recvcounts[i];

  int sendbuf_size = total_sendcodes * ( 1 /*support size*/ + dfscode_size * 5  + (dfscode_size + 1) /*max vertices*/ );
  int recvbuf_size = total_recvcodes * ( 1 /*support size*/ + dfscode_size * 5  + (dfscode_size + 1) /*max vertices*/ );

  int *sendbuf, *recvbuf;
  //int *sendbuf = new int[sendbuf_size];
  NEW_INT_ARRAY(sendbuf, sendbuf_size, *logger);
  //int *recvbuf = new int[recvbuf_size];
  NEW_INT_ARRAY(recvbuf, recvbuf_size, *logger);

  int *send_lengths = new int[numtasks];
  int *recv_lengths = new int[numtasks];

  int offset = 0;
  for(int i = 0; i<numtasks; i++) {

    send_lengths[i] = 0;

    if(rank == i)
      continue;

    for(int j = 0; j < hashed_dfscodes[i].size(); j++) {
      types::DFSCode DFS_CODE = hashed_dfscodes[i][j];
      string dfscode_str = DFS_CODE.to_string();

      sendbuf[offset++] = supports_locally_frequent[dfscode_str].size();
      send_lengths[i]++;

      for(int k = 0; k < DFS_CODE.size(); k++) {
        DFS_CODE[k].serialize( (char*) (sendbuf + offset), sizeof(int) * 5 );
        offset += 5;
        send_lengths[i] += 5;
      }

      for(int k = 0; k < supports_locally_frequent[dfscode_str].size(); k++) {
        sendbuf[offset++] = supports_locally_frequent[dfscode_str][k];
        send_lengths[i]++;
      }
    }
  }

  alltoall_msg(send_lengths, recv_lengths, 1);
  DEBUG(*logger, " send lengths = " << utils::print_array(send_lengths, numtasks) << " recv lengths = " << utils::print_array(recv_lengths, numtasks) );

  alltoallv_msg(sendbuf, recvbuf, send_lengths, recv_lengths);

  int total_send_length = offset;
  int total_recv_length = 0;
  for(int i = 0; i<numtasks; i++) {
    total_recv_length += recv_lengths[i];
  }

  int pos = 0;
  int dfscodes_received = 0;
  while(pos < total_recv_length) {
    int support_size = recvbuf[pos++];

    types::DFSCode DFS_CODE;

    for(int i = 0; i < dfscode_size; i++) {
      types::DFS dfs;
      dfs.deserialize((char*) (recvbuf + pos), 5 * sizeof(int));
      DFS_CODE.push_back(dfs);
      pos += 5;
    }

    string dfscode_str = DFS_CODE.to_string();

    DEBUG(*logger, "rank " << rank << " received hashed dfscode = " << dfscode_str);

    if(global_supports.count(dfscode_str) == 0) {
      std::vector<int> t;
      global_supports[dfscode_str] = t;

      for(int i = 0; i <support_size; i++)
        global_supports[dfscode_str].push_back(recvbuf[pos++]);

      //if could be in locally not frequent list
      if(supports_locally_not_frequent.count(dfscode_str) > 0) {
        for(int i = 0; i <support_size; i++)
          global_supports[dfscode_str][i] += supports_locally_not_frequent[dfscode_str][i];
      }

    } else {

      for(int i = 0; i <support_size; i++)
        global_supports[dfscode_str][i] += recvbuf[pos++];
    }

    dfscodes_received++;
  }

  DEBUG(*logger, " Global support size " << global_supports.size());
  delete [] sendcounts;
  delete [] recvcounts;
  DELETE_INT_ARRAY(sendbuf, *logger);
  DELETE_INT_ARRAY(recvbuf, *logger);
  //delete [] sendbuf;
  //delete [] recvbuf;
  delete [] send_lengths;
  delete [] recv_lengths;

}

void global_support_handler::hash_and_sendrecv_locally_frequent_patterns_with_int_ext_map_cardinalities(){

  // hash locally frequent dfs codes
  for(std::map<string, std::vector<int> >::iterator it = supports_locally_frequent.begin();
      it != supports_locally_frequent.end(); ++it) {

    int hashed_rank = utils::hash(it->first) % numtasks;
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);

    if(hashed_rank == rank) {            //belongs to me? initialize global support
      global_supports[it->first] = supports_locally_frequent[it->first];
      global_supports_internal[it->first] = internal_supports[it->first];
      global_supports_external[it->first] = external_supports[it->first];
      global_supports_external_max[it->first] = external_supports[it->first];
      continue;
    }

    if(hashed_dfscodes.count(hashed_rank) == 0) {
      std::vector<types::DFSCode> t;
      hashed_dfscodes[hashed_rank] = t;
    }

    hashed_dfscodes[hashed_rank].push_back(DFS_CODE);
  }

  int *sendcounts = new int[numtasks];
  int *recvcounts = new int[numtasks];
  int total_sendcodes = 0;
  int total_recvcodes = 0;

  for(int i = 0; i<numtasks; i++) {
    if(rank == i || hashed_dfscodes.count(i) == 0)
      sendcounts[i] = 0;
    else
      sendcounts[i] = hashed_dfscodes[i].size();
    total_sendcodes += sendcounts[i];
    recvcounts[i] = 0;
  }

  alltoall_msg(sendcounts, recvcounts, 1);
  DEBUG(*logger, "sent counts " << utils::print_array(sendcounts, numtasks) << " recv counts " << utils::print_array(recvcounts, numtasks));

  for(int i = 0; i<numtasks; i++)
    total_recvcodes += recvcounts[i];

  int sendbuf_size = total_sendcodes * ( 1 /*support size*/ + dfscode_size * 5  + 2 * (dfscode_size + 1) /*max vertices*/ );
  int recvbuf_size = total_recvcodes * ( 1 /*support size*/ + dfscode_size * 5  + 2 * (dfscode_size + 1) /*max vertices*/ );

  int *sendbuf, *recvbuf;
  //int *sendbuf = new int[sendbuf_size];
  NEW_INT_ARRAY(sendbuf, sendbuf_size, *logger);
  //int *recvbuf = new int[recvbuf_size];
  NEW_INT_ARRAY(recvbuf, recvbuf_size, *logger);

  int *send_lengths = new int[numtasks];
  int *recv_lengths = new int[numtasks];

  int offset = 0;
  for(int i = 0; i<numtasks; i++) {

    send_lengths[i] = 0;

    if(rank == i)
      continue;

    for(int j = 0; j < hashed_dfscodes[i].size(); j++) {
      types::DFSCode DFS_CODE = hashed_dfscodes[i][j];
      string dfscode_str = DFS_CODE.to_string();

      sendbuf[offset++] = supports_locally_frequent[dfscode_str].size();
      send_lengths[i]++;

      for(int k = 0; k < DFS_CODE.size(); k++) {
        DFS_CODE[k].serialize( (char*) (sendbuf + offset), sizeof(int) * 5 );
        offset += 5;
        send_lengths[i] += 5;
      }

      for(int k = 0; k < supports_locally_frequent[dfscode_str].size(); k++) {
        sendbuf[offset++] = internal_supports[dfscode_str][k];
        send_lengths[i]++;
        sendbuf[offset++] = external_supports[dfscode_str][k];
        send_lengths[i]++;
      }

    }
  }

  alltoall_msg(send_lengths, recv_lengths, 1);
  DEBUG(*logger, " send lengths = " << utils::print_array(send_lengths, numtasks) << " recv lengths = " << utils::print_array(recv_lengths, numtasks) );

  alltoallv_msg(sendbuf, recvbuf, send_lengths, recv_lengths);

  int total_send_length = offset;
  int total_recv_length = 0;
  for(int i = 0; i<numtasks; i++) {
    total_recv_length += recv_lengths[i];
  }

  int pos = 0;
  int dfscodes_received = 0;
  while(pos < total_recv_length) {
    int support_size = recvbuf[pos++];

    types::DFSCode DFS_CODE;

    for(int i = 0; i < dfscode_size; i++) {
      types::DFS dfs;
      dfs.deserialize((char*) (recvbuf + pos), 5 * sizeof(int));
      DFS_CODE.push_back(dfs);
      pos += 5;
    }

    string dfscode_str = DFS_CODE.to_string();

    DEBUG(*logger, "rank " << rank << " received hashed dfscode = " << dfscode_str);

    if(global_supports.count(dfscode_str) == 0) {
      std::vector<int> t;
      global_supports[dfscode_str] = t;
      global_supports_internal[dfscode_str] = t;
      global_supports_external[dfscode_str] = t;
      global_supports_external_max[dfscode_str] = t;

      for(int i = 0; i <support_size; i++) {
        global_supports_internal[dfscode_str].push_back(recvbuf[pos++]);
        global_supports_external[dfscode_str].push_back(recvbuf[pos]);
        global_supports_external_max[dfscode_str].push_back(recvbuf[pos++]);
        global_supports[dfscode_str].push_back(global_supports_internal[dfscode_str][i] + global_supports_external[dfscode_str][i]);
      }

      //if could be in locally not frequent list, in which case it was not initialized
      if(supports_locally_not_frequent.count(dfscode_str) > 0) {
        for(int i = 0; i <support_size; i++) {
          global_supports_internal[dfscode_str][i] += internal_supports[dfscode_str][i];
          global_supports_external[dfscode_str][i] += external_supports[dfscode_str][i];
          global_supports[dfscode_str][i] += internal_supports[dfscode_str][i] + external_supports[dfscode_str][i];
          if(global_supports_external_max[dfscode_str][i] < external_supports[dfscode_str][i])
            global_supports_external_max[dfscode_str][i] = external_supports[dfscode_str][i];
        }
      }

    } else {
      for(int i = 0; i <support_size; i++) {
        global_supports_internal[dfscode_str][i] += recvbuf[pos++];
        global_supports_external[dfscode_str][i] += recvbuf[pos++];
        global_supports[dfscode_str][i] += recvbuf[pos - 2] + recvbuf[pos - 1];
        if(global_supports_external_max[dfscode_str][i] < recvbuf[pos - 1])
          global_supports_external_max[dfscode_str][i] = recvbuf[pos - 1];
      }
    }

    dfscodes_received++;
  }

  DEBUG(*logger, " Global support size " << global_supports.size());

  delete [] sendcounts;
  delete [] recvcounts;

  DELETE_INT_ARRAY(sendbuf, *logger);
  DELETE_INT_ARRAY(recvbuf, *logger);
  //delete [] sendbuf;
  //delete [] recvbuf;

  delete [] send_lengths;
  delete [] recv_lengths;
}


void global_support_handler::sendrecv_unfinished_frequent_patterns(){
  int my_num_patterns = global_supports.size();
  int *all_num_patterns = new int[numtasks];
  all_gather_msg(&my_num_patterns, all_num_patterns,1);

  int *recv_lengths = new int[numtasks];
  int total_recv_len = 0;
  int total_recv_codes = 0;
  stringstream ss;
  ss << " number of patterns to consider for global support = ";
  for(int i = 0; i<numtasks; i++) {
    total_recv_codes += all_num_patterns[i];
    recv_lengths[i] = all_num_patterns[i] * dfscode_size * 5;
    ss << all_num_patterns[i] << " ";
  }

  DEBUG(*logger, ss.str());

  total_recv_len = total_recv_codes * dfscode_size * 5;
  int *send_buf, *recv_buf;
  //int *send_buf= new int[my_num_patterns * dfscode_size * 5];
  //int *recv_buf = new int[total_recv_len];
  NEW_INT_ARRAY(send_buf, my_num_patterns * dfscode_size * 5, *logger);
  NEW_INT_ARRAY(recv_buf, total_recv_len, *logger);


  //max possible for send/recv support values
  int *send_buf2, *recv_buf2;
  //int *send_buf2 = new int[total_recv_codes*(dfscode_size + 2)];
  //int *recv_buf2 = new int[total_recv_codes*(dfscode_size + 2)];
  NEW_INT_ARRAY(send_buf2, total_recv_codes * (dfscode_size + 2), *logger);
  NEW_INT_ARRAY(recv_buf2, total_recv_codes * (dfscode_size + 2), *logger);

  int *displs2 = new int[numtasks];
  int *send_lengths2 = new int[numtasks];
  int *recv_lengths2 = new int[numtasks];

  displs2[0] = 0;
  send_lengths2[0] = 0;
  recv_lengths2[0] = 0;
  for(int i = 1; i<numtasks; i++) {
    displs2[i] = displs2[i - 1] + all_num_patterns[i - 1] * (dfscode_size + 2);
    send_lengths2[i] = 0;
    recv_lengths2[i] = 0;
  }

  int offset = 0;
  //cout<<"rank " << rank <<" sending DFS codes ";
  for(std::map<string, std::vector<int> >::iterator it = global_supports.begin(); it != global_supports.end(); it++) {
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);
    for(int i = 0; i < DFS_CODE.size(); i++) {
      DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
      offset += 5;
    }
    //cout <<DFS_CODE.to_string() << " ";
  }
  //cout<<endl;

  DEBUG(*logger, " Gathering uncomputed patterns");
  //all_gatherv_msg(send_buf, recv_buf, my_num_patterns* dfscode_size * 5, recv_lengths);

  int MAX_BUFFER_SIZE = 2048;
  int max_sendrecv_count = MAX_BUFFER_SIZE / (dfscode_size * 5);
  int max_patterns_count = my_num_patterns;
  for(int i = 0; i<numtasks; i++)
    if(max_patterns_count < all_num_patterns[i])
      max_patterns_count = all_num_patterns[i];

  int iters = max_patterns_count / max_sendrecv_count;
  if(iters * max_sendrecv_count < max_patterns_count) iters++;
  int *rlens = new int[numtasks];
  int *rcounts = new int[numtasks];

  std::vector<std::vector<int> > local_supports_to_send;
  for(int i = 0; i<numtasks; i++) {
    std::vector<int> t;
    local_supports_to_send.push_back(t);
  }

  int send_offset = 0;
  for(int k = 0; k<iters; k++) {

    int send_count = max_sendrecv_count;
    if( (k + 1) * max_sendrecv_count <= my_num_patterns )
      send_count = max_sendrecv_count;
    else if (k * max_sendrecv_count < my_num_patterns )
      send_count = my_num_patterns % max_sendrecv_count;
    else
      send_count = 0;
    int send_size = send_count * dfscode_size * 5;

    for(int i = 0; i<numtasks; i++) {
      rcounts[i] = max_sendrecv_count;
      if( (k + 1) * max_sendrecv_count <= all_num_patterns[i] )
        rcounts[i] = max_sendrecv_count;
      else if (k * max_sendrecv_count < all_num_patterns[i] )
        rcounts[i] = all_num_patterns[i] % max_sendrecv_count;
      else
        rcounts[i] = 0;
      rlens[i] = rcounts[i] * dfscode_size * 5;
    }

    all_gatherv_msg(send_buf + send_offset, recv_buf, send_size, rlens);
    send_offset += send_size;

    offset = 0;
    int send_pos = 0;
    for(int i = 0; i<numtasks; i++) {
      send_pos = displs2[i] + send_lengths2[i];
      DEBUG(*logger, "RECEIVED local supports " << utils::print_array((int*) (recv_buf + offset), rlens[i]));
      if(i != rank)
        process_non_frequent_local_support_requests(recv_buf + offset, rlens[i], rcounts[i], i, send_buf2 + send_pos, send_lengths2[i]);
      else
        send_lengths2[i] += 0;
      offset += rlens[i];
    }

  }

  DEBUG(*logger, " total size = " << total_recv_codes * (dfscode_size + 2) << " displs2 = " << utils::print_array(displs2, numtasks) << " send length2 = " << utils::print_array(send_lengths2, numtasks));

  alltoall_msg(send_lengths2, recv_lengths2, 1);

  DEBUG(*logger, " send " << utils::print_array(send_lengths2, numtasks) << " recv " << utils::print_array(recv_lengths2, numtasks));


  int total_send3_length = 0;
  int total_recv3_length = 0;
  for(int i = 0; i<numtasks; i++) {
    total_send3_length += send_lengths2[i];
    total_recv3_length += recv_lengths2[i];

  }

  int *send_buf3, *recv_buf3;
  //int *send_buf3 = new int[total_send3_length];
  //int *recv_buf3 = new int[total_recv3_length];
  NEW_INT_ARRAY(send_buf3, total_send3_length, *logger);
  NEW_INT_ARRAY(recv_buf3, total_recv3_length, *logger);

  offset = 0;
  for(int i = 0; i < numtasks; i++) {
    if(i != rank)
      memcpy(send_buf3 + offset, send_buf2 + displs2[i], send_lengths2[i] * sizeof(int) );
    offset += send_lengths2[i];
  }

  alltoallv_msg(send_buf3, recv_buf3, send_lengths2, recv_lengths2);

  offset = 0;
  for(int i = 0; i<numtasks; i++) {
    if(rank != i) {
      DEBUG(*logger, "Rank " << rank << " sending to rank " << i << " " << utils::print_array(send_buf2 + displs2[i], send_lengths2[i]));
      process_received_data2(recv_buf3 + offset, recv_lengths2[i], i);
      offset += recv_lengths2[i];
    }
  }

  delete [] rlens;
  delete [] rcounts;

  delete [] all_num_patterns;
  delete [] recv_lengths;
  //delete []send_buf;
  //delete []recv_buf;
  DELETE_INT_ARRAY(send_buf, *logger);
  DELETE_INT_ARRAY(recv_buf, *logger);

  //delete []send_buf2;
  //delete []recv_buf2;
  DELETE_INT_ARRAY(send_buf2, *logger);
  DELETE_INT_ARRAY(recv_buf2, *logger);

  delete [] displs2;
  delete [] send_lengths2;
  delete [] recv_lengths2;
  //delete []send_buf3;
  //delete []recv_buf3;
  DELETE_INT_ARRAY(send_buf3, *logger);
  DELETE_INT_ARRAY(recv_buf3, *logger);


}

void global_support_handler::sendrecv_unfinished_frequent_patterns_with_int_ext_map_cardinalities(){
  int my_num_patterns = global_supports.size();
  int *all_num_patterns = new int[numtasks];
  all_gather_msg(&my_num_patterns, all_num_patterns,1);

  int *recv_lengths = new int[numtasks];
  int total_recv_len = 0;
  int total_recv_codes = 0;
  //stringstream ss;
  //ss << " number of patterns to consider for global support = ";
  for(int i = 0; i<numtasks; i++) {
    total_recv_codes += all_num_patterns[i];
    recv_lengths[i] = all_num_patterns[i] * dfscode_size * 5;
    //ss << all_num_patterns[i] << " ";
  }

  DEBUG(*logger, ss.str());

  total_recv_len = total_recv_codes * dfscode_size * 5;
  int *send_buf, *recv_buf;
  //int *send_buf= new int[my_num_patterns * dfscode_size * 5];
  //int *recv_buf = new int[total_recv_len];
  NEW_INT_ARRAY(send_buf, my_num_patterns * dfscode_size * 5, *logger);
  NEW_INT_ARRAY(recv_buf, total_recv_len, *logger);

  //max possible for send/recv support values (internal and external = times 2)
  int *send_buf2, *recv_buf2;
  //int *send_buf2 = new int[total_recv_codes*(dfscode_size + 2)];
  //int *recv_buf2 = new int[total_recv_codes*(dfscode_size + 2)];
  NEW_INT_ARRAY(send_buf2, total_recv_codes * (dfscode_size + 2) * 2, *logger);
  NEW_INT_ARRAY(recv_buf2, total_recv_codes * (dfscode_size + 2) * 2, *logger);

  int *displs2 = new int[numtasks];
  int *send_lengths2 = new int[numtasks];
  int *recv_lengths2 = new int[numtasks];

  displs2[0] = 0;
  send_lengths2[0] = 0;
  recv_lengths2[0] = 0;
  for(int i = 1; i<numtasks; i++) {
    displs2[i] = displs2[i - 1] + all_num_patterns[i - 1] * (dfscode_size + 2) * 2;
    send_lengths2[i] = 0;
    recv_lengths2[i] = 0;
  }

  int offset = 0;
  //cout<<"rank " << rank <<" sending DFS codes ";
  for(std::map<string, std::vector<int> >::iterator it = global_supports.begin(); it != global_supports.end(); it++) {
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);
    for(int i = 0; i < DFS_CODE.size(); i++) {
      DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
      offset += 5;
    }
    //cout <<DFS_CODE.to_string() << " ";
  }
  //cout<<endl;

  DEBUG(*logger, " Gathering uncomputed patterns");
  //all_gatherv_msg(send_buf, recv_buf, my_num_patterns* dfscode_size * 5, recv_lengths);

  int MAX_BUFFER_SIZE = 2048;
  int max_sendrecv_count = MAX_BUFFER_SIZE / (dfscode_size * 5);
  int max_patterns_count = my_num_patterns;
  for(int i = 0; i<numtasks; i++)
    if(max_patterns_count < all_num_patterns[i])
      max_patterns_count = all_num_patterns[i];

  int iters = max_patterns_count / max_sendrecv_count;
  if(iters * max_sendrecv_count < max_patterns_count) iters++;
  int *rlens = new int[numtasks];
  int *rcounts = new int[numtasks];

  std::vector<std::vector<int> > local_supports_to_send;
  for(int i = 0; i<numtasks; i++) {
    std::vector<int> t;
    local_supports_to_send.push_back(t);
  }

  int send_offset = 0;
  for(int k = 0; k<iters; k++) {

    int send_count = max_sendrecv_count;
    if( (k + 1) * max_sendrecv_count <= my_num_patterns )
      send_count = max_sendrecv_count;
    else if (k * max_sendrecv_count < my_num_patterns )
      send_count = my_num_patterns % max_sendrecv_count;
    else
      send_count = 0;
    int send_size = send_count * dfscode_size * 5;

    for(int i = 0; i<numtasks; i++) {
      rcounts[i] = max_sendrecv_count;
      if( (k + 1) * max_sendrecv_count <= all_num_patterns[i] )
        rcounts[i] = max_sendrecv_count;
      else if (k * max_sendrecv_count < all_num_patterns[i] )
        rcounts[i] = all_num_patterns[i] % max_sendrecv_count;
      else
        rcounts[i] = 0;
      rlens[i] = rcounts[i] * dfscode_size * 5;
    }

    all_gatherv_msg(send_buf + send_offset, recv_buf, send_size, rlens);
    send_offset += send_size;

    offset = 0;
    int send_pos = 0;
    for(int i = 0; i<numtasks; i++) {
      send_pos = displs2[i] + send_lengths2[i];
      DEBUG(*logger, "RECEIVED local supports " << utils::print_array((int*) (recv_buf + offset), rlens[i]));
      if(i != rank)
        process_non_frequent_local_support_requests_with_int_ext_map_cardinalities(recv_buf + offset, rlens[i], rcounts[i], i, send_buf2 + send_pos, send_lengths2[i]);
      else
        send_lengths2[i] += 0;

      offset += rlens[i];
    }
  }

  DEBUG(*logger, " total size = " << total_recv_codes * (dfscode_size + 2) << " displs2 = " << utils::print_array(displs2, numtasks) << " send length2 = " << utils::print_array(send_lengths2, numtasks));

  alltoall_msg(send_lengths2, recv_lengths2, 1);

  DEBUG(*logger, " send " << utils::print_array(send_lengths2, numtasks) << " recv " << utils::print_array(recv_lengths2, numtasks));

  int total_send3_length = 0;
  int total_recv3_length = 0;
  for(int i = 0; i<numtasks; i++) {
    total_send3_length += send_lengths2[i];
    total_recv3_length += recv_lengths2[i];
  }

  int *send_buf3, *recv_buf3;
  //int *send_buf3 = new int[total_send3_length];
  //int *recv_buf3 = new int[total_recv3_length];
  NEW_INT_ARRAY(send_buf3, total_send3_length, *logger);
  NEW_INT_ARRAY(recv_buf3, total_recv3_length, *logger);

  offset = 0;
  for(int i = 0; i < numtasks; i++) {
    if(i != rank)
      memcpy(send_buf3 + offset, send_buf2 + displs2[i], send_lengths2[i] * sizeof(int) );
    offset += send_lengths2[i];
  }

  alltoallv_msg(send_buf3, recv_buf3, send_lengths2, recv_lengths2);

  offset = 0;
  for(int i = 0; i<numtasks; i++) {
    if(rank != i) {
      DEBUG(*logger, "Rank " << rank << " sending to rank " << i << " " << utils::print_array(send_buf2 + displs2[i], send_lengths2[i]));
      process_received_data2_with_int_ext_map_cardinalities(recv_buf3 + offset, recv_lengths2[i], i);
      offset += recv_lengths2[i];
    }
  }

  delete [] rlens;
  delete [] rcounts;

  delete [] all_num_patterns;
  delete [] recv_lengths;
  //delete []send_buf;
  //delete []recv_buf;
  DELETE_INT_ARRAY(send_buf, *logger);
  DELETE_INT_ARRAY(recv_buf, *logger);

  //delete []send_buf2;
  //delete []recv_buf2;
  DELETE_INT_ARRAY(send_buf2, *logger);
  DELETE_INT_ARRAY(recv_buf2, *logger);

  delete [] displs2;
  delete [] send_lengths2;
  delete [] recv_lengths2;
  //delete []send_buf3;
  //delete []recv_buf3;
  DELETE_INT_ARRAY(send_buf3, *logger);
  DELETE_INT_ARRAY(recv_buf3, *logger);

}



void global_support_handler::sendrecv_unfinished_frequent_patterns_with_int_ext_map_cardinalities2(){

  //Every process their remaining patterns, i.e., not determined to be frequent
  int my_num_patterns = global_supports.size();
  int *all_num_patterns = new int[numtasks];
  all_gather_msg(&my_num_patterns, all_num_patterns,1);

  int *recv_lengths = new int[numtasks];
  int total_recv_len = 0;
  int total_recv_codes = 0;
  stringstream ss;
  ss << " number of patterns to consider for global support = ";
  for(int i = 0; i<numtasks; i++) {
    total_recv_codes += all_num_patterns[i];
    recv_lengths[i] = all_num_patterns[i] * dfscode_size * 5;
    ss << all_num_patterns[i] << " ";
  }

  INFO(*logger, ss.str());

  total_recv_len = total_recv_codes * dfscode_size * 5;
  //int *send_buf, *recv_buf;
  int *send_buf = new int[my_num_patterns * dfscode_size * 5];
  int *recv_buf = new int[total_recv_len];
  //NEW_INT_ARRAY(send_buf, my_num_patterns * dfscode_size * 5, *logger);
  //NEW_INT_ARRAY(recv_buf, total_recv_len, *logger);

  //max possible for send/recv support values (internal and external = times 2)
  //int *send_buf2, *recv_buf2;
  int *send_buf2 = new int[total_recv_codes * (dfscode_size + 2) * 2];
  int *recv_buf2 = new int[total_recv_codes * (dfscode_size + 2) * 2];
  //NEW_INT_ARRAY(send_buf2, total_recv_codes * (dfscode_size + 2) * 2, *logger);
  //NEW_INT_ARRAY(recv_buf2, total_recv_codes * (dfscode_size + 2) * 2, *logger);

  int *displs2 = new int[numtasks];
  int *send_lengths2 = new int[numtasks];
  int *recv_lengths2 = new int[numtasks];

  displs2[0] = 0;
  send_lengths2[0] = 0;
  recv_lengths2[0] = 0;
  for(int i = 1; i<numtasks; i++) {
    displs2[i] = displs2[i - 1] + all_num_patterns[i - 1] * (dfscode_size + 2) * 2;
    send_lengths2[i] = 0;
    recv_lengths2[i] = 0;
  }

  int offset = 0;
  //cout<<"rank " << rank <<" sending DFS codes ";
  /*for(std::map<string, std::vector<int> >::iterator it = global_supports.begin(); it != global_supports.end(); it++) {
     types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);
     for(int i = 0; i < DFS_CODE.size(); i++) {
      DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
      offset += 5;
     }
     //cout <<DFS_CODE.to_string() << " ";
     }*/
  //cout<<endl;

  //Parallel serialization (from above)
  int numth = numthreads;
  if(numth > global_supports.size())
    numth = global_supports.size();

  cout << "NUMBER OF threads = " << numth << endl;
#pragma omp parallel num_threads(numth) private(offset)
  {
    int thread_id = omp_get_thread_num();
    int dfscodes_per_thread = (int) ceil(global_supports.size() * 1.0 / numth);
    int start_index = thread_id * dfscodes_per_thread;
    int end_index = start_index + dfscodes_per_thread - 1;

    if (end_index > global_supports.size() - 1)
      end_index = global_supports.size() - 1;

    //omp_set_lock(&lock_gs);
    //cout<<"rank "<<rank<<"thread"<<thread_id<<" start = "<<start_index << " end = " << end_index <<endl;
    //omp_unset_lock(&lock_gs);

    offset = start_index * 5 * dfscode_size;
    std::map<string, std::vector<int> >::iterator it = global_supports.begin();
    std::advance( it, start_index );
    for(int k = start_index; k <= end_index; k++, it++) {
      //for(std::map<string, std::vector<int> >::iterator it = global_supports.begin() + start_index ; it != global_supports.begin() + end_index + 1; it++) {
      types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);
      for(int i = 0; i < DFS_CODE.size(); i++) {
        DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
        offset += 5;
      }
    }
  } //pragma



  DEBUG(*logger, " Gathering uncomputed patterns");
  //all_gatherv_msg(send_buf, recv_buf, my_num_patterns* dfscode_size * 5, recv_lengths);

  int MAX_BUFFER_SIZE = 20480 * 5;
  int max_sendrecv_count = MAX_BUFFER_SIZE / (dfscode_size * 5);
  int max_patterns_count = my_num_patterns;
  for(int i = 0; i<numtasks; i++)
    if(max_patterns_count < all_num_patterns[i])
      max_patterns_count = all_num_patterns[i];

  int iters = max_patterns_count / max_sendrecv_count;
  if(iters * max_sendrecv_count < max_patterns_count) iters++;
  int *rlens = new int[numtasks];
  int *rcounts = new int[numtasks];

  std::vector<std::vector<int> > local_supports_to_send;
  for(int i = 0; i<numtasks; i++) {
    std::vector<int> t;
    local_supports_to_send.push_back(t);
  }

  int send_offset = 0;
  for(int k = 0; k<iters; k++) {

    int send_count = max_sendrecv_count;
    if( (k + 1) * max_sendrecv_count <= my_num_patterns )
      send_count = max_sendrecv_count;
    else if (k * max_sendrecv_count < my_num_patterns )
      send_count = my_num_patterns % max_sendrecv_count;
    else
      send_count = 0;
    int send_size = send_count * dfscode_size * 5;

    for(int i = 0; i<numtasks; i++) {
      rcounts[i] = max_sendrecv_count;
      if( (k + 1) * max_sendrecv_count <= all_num_patterns[i] )
        rcounts[i] = max_sendrecv_count;
      else if (k * max_sendrecv_count < all_num_patterns[i] )
        rcounts[i] = all_num_patterns[i] % max_sendrecv_count;
      else
        rcounts[i] = 0;
      rlens[i] = rcounts[i] * dfscode_size * 5;
    }

    all_gatherv_msg(send_buf + send_offset, recv_buf, send_size, rlens);
    send_offset += send_size;

    offset = 0;
    int send_pos = 0;
    //serial version
    /*for(int i = 0; i<numtasks; i++) {
       send_pos = displs2[i] + send_lengths2[i];
       DEBUG(*logger, "RECEIVED local supports " << utils::print_array((int*) (recv_buf + offset), rlens[i]));
       if(i != rank)
        process_non_frequent_local_support_requests_with_int_ext_map_cardinalities(recv_buf + offset, rlens[i], rcounts[i], i, send_buf2 + send_pos, send_lengths2[i]);
       else
        send_lengths2[i] += 0;
       offset += rlens[i];
       }*/

    //parallelize above code
    numth = numthreads;
    if(numth > numtasks)
      numth = numtasks;     //no parallelism

    cout << "NUMBER OF threads = " << numth << endl;
#pragma omp parallel for num_threads(numth) private(offset)
    for(int i = 0; i<numtasks; i++) {

      offset = 0;
      if(i > 0 )
        for(int j = 0; j < i; j++)
          offset += rlens[j];

      send_pos = displs2[i] + send_lengths2[i];
      DEBUG(*logger, "RECEIVED local supports " << utils::print_array((int*) (recv_buf + offset), rlens[i]));
      if(i != rank)
        process_non_frequent_local_support_requests_with_int_ext_map_cardinalities(recv_buf + offset, rlens[i], rcounts[i], i, send_buf2 + send_pos, send_lengths2[i]);
      else
        send_lengths2[i] += 0;
    }


  }

  DEBUG(*logger, " total size = " << total_recv_codes * (dfscode_size + 2) << " displs2 = " << utils::print_array(displs2, numtasks) << " send length2 = " << utils::print_array(send_lengths2, numtasks));

  alltoall_msg(send_lengths2, recv_lengths2, 1);

  DEBUG(*logger, " send " << utils::print_array(send_lengths2, numtasks) << " recv " << utils::print_array(recv_lengths2, numtasks));

  int total_send3_length = 0;
  int total_recv3_length = 0;
  for(int i = 0; i<numtasks; i++) {
    total_send3_length += send_lengths2[i];
    total_recv3_length += recv_lengths2[i];
  }

  //int *send_buf3, *recv_buf3;
  int *send_buf3 = new int[total_send3_length];
  int *recv_buf3 = new int[total_recv3_length];
  //NEW_INT_ARRAY(send_buf3, total_send3_length, *logger);
  //NEW_INT_ARRAY(recv_buf3, total_recv3_length, *logger);

  offset = 0;
  for(int i = 0; i < numtasks; i++) {
    if(i != rank)
      memcpy(send_buf3 + offset, send_buf2 + displs2[i], send_lengths2[i] * sizeof(int) );
    offset += send_lengths2[i];
  }

  //alltoallv_msg(send_buf2, recv_buf2, send_lengths2, recv_lengths2, displs2, displs2);

  alltoallv_msg(send_buf3, recv_buf3, send_lengths2, recv_lengths2);

  //serial version
  /*offset = 0;
     for(int i = 0; i<numtasks; i++) {
     if(rank != i) {
      DEBUG(*logger, "Rank " << rank << " sending to rank " << i << " " << utils::print_array(send_buf2 + displs2[i], send_lengths2[i]));
      //process_received_data2(recv_buf2 + displs2[i], recv_lengths2[i], i);
      process_received_data2_with_int_ext_map_cardinalities(recv_buf3 + offset, recv_lengths2[i], i);
      offset += recv_lengths2[i];
     }
     }*/

  //parallelize above code
  numth = numthreads;
  if(numth > numtasks)
    numth = numtasks;       //no parallelism
  cout << "NUMBER OF threads = " << numth << endl;
#pragma omp parallel for num_threads(numth) private(offset)
  for(int i = 0; i<numtasks; i++) {
    offset = 0;
    if(i > 0 )
      for(int j = 0; j < i; j++)
        offset += recv_lengths2[j];

    if(rank != i) {
      DEBUG(*logger, "Rank " << rank << " sending to rank " << i << " " << utils::print_array(send_buf2 + displs2[i], send_lengths2[i]));
      process_received_data2_with_int_ext_map_cardinalities(recv_buf3 + offset, recv_lengths2[i], i);
    }
  }


  delete [] rlens;
  delete [] rcounts;

  delete [] all_num_patterns;
  delete [] recv_lengths;
  delete [] send_buf;
  delete [] recv_buf;
  //DELETE_INT_ARRAY(send_buf, *logger);
  //DELETE_INT_ARRAY(recv_buf, *logger);

  delete [] send_buf2;
  delete [] recv_buf2;
  //DELETE_INT_ARRAY(send_buf2, *logger);
  //DELETE_INT_ARRAY(recv_buf2, *logger);

  delete [] displs2;
  delete [] send_lengths2;
  delete [] recv_lengths2;
  delete [] send_buf3;
  delete [] recv_buf3;
  //DELETE_INT_ARRAY(send_buf3, *logger);
  //DELETE_INT_ARRAY(recv_buf3, *logger);

}

void global_support_handler::get_all_globally_estimated_patterns(std::vector<types::DFSCode> &globally_frequent){
  globally_frequent = all_globally_estimated_dfscodes;
}

void global_support_handler::remove_all_globally_estimated_patterns(){
  all_globally_estimated_dfscodes.clear();
}

void global_support_handler::get_all_globally_estimated_patterns(std::vector<types::DFSCode> &globally_frequent, int max_support){
  for(int i = 0; i < all_globally_estimated_dfscodes.size(); i++) {
    if(final_dfscode_supports[all_globally_estimated_dfscodes[i].to_string()] <= max_support)
      globally_frequent.push_back(all_globally_estimated_dfscodes[i]);
  }
}

bool global_support_handler::remove_globally_frequent_dfscode(types::DFSCode dfscode){
  //makes sure the count decremented only once
  int hashed_rank = utils::hash(dfscode.to_string()) % numtasks;
  if(rank == hashed_rank)
    frequent_patterns--;

  for(std::vector<types::DFSCode>::iterator it = globally_frequent_dfscodes.begin();
      it != globally_frequent_dfscodes.end(); it++) {
    if(*it == dfscode) {
      DEBUG(*logger, "DFSCODE found " << dfscode);
      globally_frequent_dfscodes.erase(it);
      //hashed_rank = hashed_rank % num_partitions;
      return true;
    }
  }
  return false;
}

void global_support_handler::set_dfscode_support(types::DFSCode dfscode, int support){
  final_dfscode_supports[dfscode.to_string()] = support;
}

void global_support_handler::send_globally_frequent_patterns(){
  int my_num_patterns = globally_frequent_dfscodes_hashed.size();
  int *all_num_patterns = new int[numtasks];
  all_gather_msg(&my_num_patterns, all_num_patterns,1);

  int *recv_lengths = new int[numtasks];
  int total_recv_len = 0;
  stringstream ss;
  ss << " number of globally frequent patterns = ";
  for(int i = 0; i<numtasks; i++) {
    total_recv_len += all_num_patterns[i] * dfscode_size * 5;
    recv_lengths[i] = all_num_patterns[i] * dfscode_size * 5;
    ss << all_num_patterns[i] << " ";
  }

  DEBUG(*logger, ss.str());

  int *send_buf = new int[my_num_patterns * dfscode_size * 5];
  int *recv_buf = new int[total_recv_len];

  int offset = 0;
  //cout<<"rank " << rank <<" sending DFS codes ";
  for(int k = 0; k < globally_frequent_dfscodes_hashed.size(); k++) {
    types::DFSCode DFS_CODE = globally_frequent_dfscodes_hashed[k];
    for(int i = 0; i < DFS_CODE.size(); i++) {
      DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
      offset += 5;
    }
    //cout <<DFS_CODE.to_string() << " ";

    //insert globally frequent pattern into my list
    if(supports_locally_frequent.count(DFS_CODE.to_string()) > 0 ||
       supports_locally_not_frequent.count(DFS_CODE.to_string()) > 0) {
      globally_frequent_dfscodes.push_back(DFS_CODE);
    }
    all_globally_estimated_dfscodes.push_back(DFS_CODE);
  }
  //cout<<endl;

  all_gatherv_msg(send_buf, recv_buf, my_num_patterns * dfscode_size * 5, recv_lengths);

  ss.clear();
  ss << "rank " << rank << " Received globally frequent dfscodes ";
  offset = 0;
  for(int k = 0; k<numtasks; k++) {
    for(int i = 0; i<recv_lengths[k]; i++)
      ss << recv_buf[offset + i] <<  " ";
    ss << " | ";
    offset += recv_lengths[k];
  }
  DEBUG(*logger, ss.str());

  offset = 0;
  for(int i = 0; i<numtasks; i++) {
    if(i != rank)
      process_globally_frequent_patterns(recv_buf + offset, recv_lengths[i], all_num_patterns[i], i);
    offset += recv_lengths[i];
  }

  delete [] all_num_patterns;
  delete [] send_buf;
  delete [] recv_buf;
  delete [] recv_lengths;
}


void global_support_handler::send_globally_frequent_and_estimated_patterns(){
  int my_num_patterns = globally_frequent_dfscodes_hashed.size();
  int *all_num_patterns = new int[numtasks];
  all_gather_msg(&my_num_patterns, all_num_patterns,1);

  int *recv_lengths = new int[numtasks];
  int total_recv_len = 0;
  stringstream ss;
  ss << " number of globally frequent patterns = ";
  for(int i = 0; i<numtasks; i++) {
    total_recv_len += all_num_patterns[i] * dfscode_size * 5 + all_num_patterns[i];            // second one is the flag to indicate whether it is actually globally frequent or not
    recv_lengths[i] = all_num_patterns[i] * dfscode_size * 5 + all_num_patterns[i];
    ss << all_num_patterns[i] << " ";
  }

  DEBUG(*logger, ss.str());

  int *send_buf = new int[my_num_patterns * dfscode_size * 5 + my_num_patterns];
  int *recv_buf = new int[total_recv_len];

  int offset = 0;
  int marked_freq = 0;
  //cout<<"rank " << rank <<" sending DFS codes ";
  for(int k = 0; k < globally_frequent_dfscodes_hashed.size(); k++) {
    types::DFSCode DFS_CODE = globally_frequent_dfscodes_hashed[k];
    for(int i = 0; i < DFS_CODE.size(); i++) {
      DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
      offset += 5;
    }
    //cout <<DFS_CODE.to_string() << " ";
    send_buf[offset++] = marked_actual_globally_frequent_patterns[DFS_CODE.to_string()];
    if( marked_actual_globally_frequent_patterns[DFS_CODE.to_string()]  == 1)
      marked_freq++;

    //insert globally frequent pattern into my list
    if(supports_locally_frequent.count(DFS_CODE.to_string()) > 0 || supports_locally_not_frequent.count(DFS_CODE.to_string()) > 0) {
      globally_frequent_dfscodes.push_back(DFS_CODE);
    }
    if(marked_actual_globally_frequent_patterns[DFS_CODE.to_string()] == 0)
      all_globally_estimated_dfscodes.push_back(DFS_CODE);
  }
  //cout<<endl;
  DEBUG(*logger, "Total estimated frequent = " << marked_actual_globally_frequent_patterns.size() << " actually freq = " << marked_freq);

  all_gatherv_msg(send_buf, recv_buf, my_num_patterns * dfscode_size * 5 + my_num_patterns, recv_lengths);

  ss.clear();
  ss << "rank " << rank << " Received globally frequent dfscodes ";
  offset = 0;
  for(int k = 0; k<numtasks; k++) {
    for(int i = 0; i<recv_lengths[k]; i++)
      ss << recv_buf[offset + i] <<  " ";
    ss << " | ";
    offset += recv_lengths[k];
  }
  DEBUG(*logger, ss.str());

  offset = 0;
  for(int i = 0; i<numtasks; i++) {
    if(i != rank)
      process_globally_frequent_patterns_with_frequent_flags(recv_buf + offset, recv_lengths[i], all_num_patterns[i], i);
    offset += recv_lengths[i];
  }

  delete [] all_num_patterns;
  delete [] send_buf;
  delete [] recv_buf;
  delete [] recv_lengths;
}

void global_support_handler::send_globally_frequent_patterns_to_rank_group(){

  MPI_Comm rank_group;
  int membershipKey = rank % num_partitions;
  MPI_Comm_split(MPI_COMM_WORLD, membershipKey, rank, &rank_group);

  int group_size = numtasks / num_partitions;
  int my_num_patterns = globally_frequent_dfscodes.size();
  int *all_num_patterns = new int[group_size];
  all_gather_msg(&my_num_patterns, all_num_patterns, 1, rank_group);

  int *recv_lengths = new int[group_size];
  int total_recv_len = 0;
  stringstream ss;
  ss << "From rank group: number of globally frequent patterns = ";
  for(int i = 0; i<group_size; i++) {
    total_recv_len += all_num_patterns[i] * dfscode_size * 5;
    recv_lengths[i] = all_num_patterns[i] * dfscode_size * 5;
    ss << all_num_patterns[i] << " ";
  }

  DEBUG(*logger, ss.str());

  int *send_buf, *recv_buf;
  //int *send_buf = new int[my_num_patterns * dfscode_size * 5];
  //int *recv_buf = new int[total_recv_len];
  NEW_INT_ARRAY(send_buf, my_num_patterns * dfscode_size * 5, *logger);
  NEW_INT_ARRAY(recv_buf, total_recv_len, *logger);

  int offset = 0;
  //cout<<"rank " << rank <<" sending DFS codes ";
  for(int k = 0; k < globally_frequent_dfscodes.size(); k++) {
    types::DFSCode DFS_CODE = globally_frequent_dfscodes[k];
    for(int i = 0; i < DFS_CODE.size(); i++) {
      DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
      offset += 5;
    }
    //cout <<DFS_CODE.to_string() << " ";
  }
  //cout<<endl;

  all_gatherv_msg(send_buf, recv_buf, my_num_patterns * dfscode_size * 5, recv_lengths, rank_group);


  //remove current list of dfscods so all ranks in the same group have the same order of the codes
  //will be refilled automatically
  globally_frequent_dfscodes.clear();

  offset = 0;
  for(int i = 0; i<group_size; i++) {
    process_globally_frequent_patterns2(recv_buf + offset, recv_lengths[i], all_num_patterns[i]);
    offset += recv_lengths[i];
  }

  delete [] all_num_patterns;
  //delete []send_buf;
  //delete []recv_buf;
  DELETE_INT_ARRAY(send_buf, *logger);
  DELETE_INT_ARRAY(recv_buf, *logger);

  delete [] recv_lengths;
}

void global_support_handler::compute_globally_frequent_patterns(){

  DEBUG(*logger, "Global support dfscode count = " << global_supports.size());
  for(std::map<string, std::vector<int> >::iterator it = global_supports.begin(); it != global_supports.end(); ) {

    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);

    DEBUG(*logger, "global support of " << it->first << " = " << utils::print_vector(global_supports[it->first]));
    unsigned int global_support = 0xFFFFFFFF;

    for(int i = 0; i< (it->second).size(); i++) {
      if(global_support > (it->second)[i])
        global_support = (it->second)[i];
    }

    if(global_support >= minsupp) {
      //report DFS_CODE
      frequent_patterns++;
      globally_frequent_dfscodes_hashed.push_back(DFS_CODE);
      final_dfscode_supports[DFS_CODE.to_string()] = global_support;
      /*if(supports_locally_frequent.count(DFS_CODE.to_string()) > 0)
              local_supports_globally_frequent[DFS_CODE.to_string()] = *std::min_element(supports_locally_frequent[DFS_CODE.to_string()].begin(), supports_locally_frequent[DFS_CODE.to_string()].end());
         else
              local_supports_globally_frequent[DFS_CODE.to_string()] = *std::min_element(supports_locally_not_frequent[DFS_CODE.to_string()].begin(), supports_locally_not_frequent[DFS_CODE.to_string()].end());
       */
      //if(out)
      //out->output_graph(DFS_CODE, global_support);
      global_supports.erase(it++);
    }else{
      ++it;
    }
  }
  DEBUG(*logger, "Global support dfscode count = " << global_supports.size() << " frequent " << frequent_patterns);
}


void global_support_handler::compute_globally_frequent_patterns_with_int_ext_map_cardinalities(){

  DEBUG(*logger, "Global support dfscode count = " << global_supports.size());
  int original = global_supports.size(), actual_frequent = 0, candidates = 0;

  for(std::map<string, std::vector<int> >::iterator it = global_supports.begin(); it != global_supports.end(); ) {

    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);

    DEBUG(*logger, "global support of " << it->first << " = " << utils::print_vector(global_supports[it->first]));

    unsigned int global_support_internal_external = 0xFFFFFFFF;
    for(int i = 0; i< (it->second).size(); i++) {
      int max_int_ext = global_supports_internal[it->first][i];
      if(global_support_internal_external > max_int_ext)
        global_support_internal_external = max_int_ext;
    }

    unsigned int global_support = 0xFFFFFFFF;

    for(int i = 0; i< (it->second).size(); i++) {
      if(global_support > (it->second)[i] )
        global_support = (it->second)[i];
    }

    if(global_support_internal_external >= minsupp) {
      //report DFS_CODE
      frequent_patterns++;
      actual_frequent++;
      globally_frequent_dfscodes_hashed.push_back(DFS_CODE);
      marked_actual_globally_frequent_patterns[DFS_CODE.to_string()] = 1;
      final_dfscode_supports[DFS_CODE.to_string()] = global_support;
      global_supports.erase(it++);
    }else{
      marked_actual_globally_frequent_patterns[DFS_CODE.to_string()] = 0;
      if(global_support >= minsupp) {
        //report DFS_CODE
        frequent_patterns++;
        candidates++;
        globally_frequent_dfscodes_hashed.push_back(DFS_CODE);
        final_dfscode_supports[DFS_CODE.to_string()] = global_support;
        global_supports.erase(it++);
      }else{
        ++it;
      }
    }
  }

  DEBUG(*logger, "Global support dfscode count = " << global_supports.size() << " frequent " << frequent_patterns);
  DEBUG(*logger, "rank = " << rank << " level = " << dfscode_size << " all patterns = " << global_supports.size() << " actually freq = " << actual_frequent << " other estimated = " << candidates << " total = " << (actual_frequent + candidates) << " cumulative frequent (incl false pos for this level) = " << frequent_patterns );
}


void global_support_handler::get_all_locally_frequent_patterns(std::vector<types::DFSCode> &all_locally_frequent){

  std::set<string> locally_frequent_set;

  int my_num_patterns = supports_locally_frequent.size();
  int *all_num_patterns = new int[numtasks];
  all_gather_msg(&my_num_patterns, all_num_patterns, 1);

  int *recv_lengths = new int[numtasks];
  int total_recv_len = 0;
  stringstream ss;
  ss << " number of globally frequent patterns = ";
  for(int i = 0; i<numtasks; i++) {
    total_recv_len += all_num_patterns[i] * dfscode_size * 5;
    recv_lengths[i] = all_num_patterns[i] * dfscode_size * 5;
    ss << all_num_patterns[i] << " ";
  }

  DEBUG(*logger, ss.str());

  int *send_buf = new int[my_num_patterns * dfscode_size * 5];
  int *recv_buf = new int[total_recv_len];

  int offset = 0;
  //cout<<"rank " << rank <<" sending DFS codes ";
  for(std::map<string, std::vector<int> >::iterator it = supports_locally_frequent.begin(); it != supports_locally_frequent.end(); it++) {
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);
    for(int i = 0; i < DFS_CODE.size(); i++) {
      DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
      offset += 5;
    }
    //cout <<DFS_CODE.to_string() << " ";

    locally_frequent_set.insert(it->first);
    //all_locally_frequent.push_back(DFS_CODE);
  }
  //cout<<endl;

  all_gatherv_msg(send_buf, recv_buf, my_num_patterns * dfscode_size * 5, recv_lengths);

  ss.clear();
  ss << "rank " << rank << " Received locally frequent dfscodes ";
  offset = 0;
  for(int k = 0; k<numtasks; k++) {
    for(int i = 0; i<recv_lengths[k]; i++)
      ss << recv_buf[offset + i] <<  " ";
    ss << " | ";
    offset += recv_lengths[k];
  }
  DEBUG(*logger, ss.str());

  offset = 0;
  for(int i = 0; i<numtasks; i++) {
    if(i != rank) {
      //process_globally_frequent_patterns(recv_buf + offset, recv_lengths[i], all_num_patterns[i], i);
      int *buffer = recv_buf + offset;
      int len = recv_lengths[i];
      int recv_pos = 0;
      while(recv_pos < len) {
        //get DFS code
        types::DFSCode DFS_CODE;
        for(int j = 0; j < dfscode_size; j++) {
          types::DFS dfs;
          dfs.deserialize((char*) (buffer + recv_pos), 5 * sizeof(int));
          DFS_CODE.push_back(dfs);
          recv_pos += 5;
        }

        DEBUG(*logger, "Locally frequent dfscode received " << DFS_CODE.to_string());
        locally_frequent_set.insert(DFS_CODE.to_string());
      }

    }
    offset += recv_lengths[i];
  }


  for(std::set<string>::iterator it = locally_frequent_set.begin(); it != locally_frequent_set.end(); it++) {
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(*it);
    all_locally_frequent.push_back(DFS_CODE);
  }

  DEBUG(*logger, "Total locally frequent = " << all_locally_frequent.size());
  delete [] all_num_patterns;
  delete [] send_buf;
  delete [] recv_buf;
  delete [] recv_lengths;

}



void global_support_handler::get_all_locally_not_frequent_patterns(std::vector<types::DFSCode> &all_locally_not_frequent){

  std::set<string> locally_not_frequent_set;

  int my_num_patterns = supports_locally_not_frequent.size();
  int *all_num_patterns = new int[numtasks];
  all_gather_msg(&my_num_patterns, all_num_patterns, 1);

  int *recv_lengths = new int[numtasks];
  int total_recv_len = 0;
  stringstream ss;
  ss << " number of globally frequent patterns = ";
  for(int i = 0; i<numtasks; i++) {
    total_recv_len += all_num_patterns[i] * dfscode_size * 5;
    recv_lengths[i] = all_num_patterns[i] * dfscode_size * 5;
    ss << all_num_patterns[i] << " ";
  }

  DEBUG(*logger, ss.str());

  int *send_buf = new int[my_num_patterns * dfscode_size * 5];
  int *recv_buf = new int[total_recv_len];

  int offset = 0;
  //cout<<"rank " << rank <<" sending DFS codes ";
  for(std::map<string, std::vector<int> >::iterator it = supports_locally_not_frequent.begin(); it != supports_locally_not_frequent.end(); it++) {
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);
    for(int i = 0; i < DFS_CODE.size(); i++) {
      DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
      offset += 5;
    }
    //cout <<DFS_CODE.to_string() << " ";

    locally_not_frequent_set.insert(it->first);
    //all_locally_frequent.push_back(DFS_CODE);
  }
  //cout<<endl;

  all_gatherv_msg(send_buf, recv_buf, my_num_patterns * dfscode_size * 5, recv_lengths);

  ss.clear();
  ss << "rank " << rank << " Received locally not frequent dfscodes ";
  offset = 0;
  for(int k = 0; k<numtasks; k++) {
    for(int i = 0; i<recv_lengths[k]; i++)
      ss << recv_buf[offset + i] <<  " ";
    ss << " | ";
    offset += recv_lengths[k];
  }
  DEBUG(*logger, ss.str());

  offset = 0;
  for(int i = 0; i<numtasks; i++) {
    if(i != rank) {
      //process_globally_frequent_patterns(recv_buf + offset, recv_lengths[i], all_num_patterns[i], i);
      int *buffer = recv_buf + offset;
      int len = recv_lengths[i];
      int recv_pos = 0;
      while(recv_pos < len) {
        //get DFS code
        types::DFSCode DFS_CODE;
        for(int j = 0; j < dfscode_size; j++) {
          types::DFS dfs;
          dfs.deserialize((char*) (buffer + recv_pos), 5 * sizeof(int));
          DFS_CODE.push_back(dfs);
          recv_pos += 5;
        }

        DEBUG(*logger, "Locally not frequent dfscode received " << DFS_CODE.to_string());
        locally_not_frequent_set.insert(DFS_CODE.to_string());
      }

    }
    offset += recv_lengths[i];
  }


  std::vector<types::DFSCode> all_locally_frequent;
  get_all_locally_frequent_patterns(all_locally_frequent);
  std::set<string> all_locally_frequent_set;
  for(std::vector<types::DFSCode>::iterator it = all_locally_frequent.begin(); it != all_locally_frequent.end(); it++)
    all_locally_frequent_set.insert(it->to_string());

  for(std::set<string>::iterator it = locally_not_frequent_set.begin(); it != locally_not_frequent_set.end(); it++) {
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(*it);

    if(all_locally_frequent_set.find(*it) == all_locally_frequent_set.end())
      all_locally_not_frequent.push_back(DFS_CODE);
  }

  if(rank == 0) {
    INFO(*logger, "level = " << dfscode_size << " Total locally frequent = " << all_locally_frequent.size());
    INFO(*logger, "level = " << dfscode_size << " Total locally NOT frequent = " << all_locally_not_frequent.size() << " SET size was " << locally_not_frequent_set.size());
  }
  delete [] all_num_patterns;
  delete [] send_buf;
  delete [] recv_buf;
  delete [] recv_lengths;
}



void global_support_handler::gather_locally_frequent_patterns_with_no_estimation(int dfscode_size){

  this->dfscode_size = dfscode_size;
  get_all_locally_frequent_patterns(all_globally_estimated_dfscodes);

  if(rank == 0)
    frequent_patterns += all_globally_estimated_dfscodes.size();

  for(std::vector<types::DFSCode>::iterator it = all_globally_estimated_dfscodes.begin(); it != all_globally_estimated_dfscodes.end(); it++) {
    if(supports_locally_frequent.count(it->to_string()) > 0 || supports_locally_not_frequent.count(it->to_string()) > 0) {
      globally_frequent_dfscodes.push_back(*it);
    }
  }

  INFO(*logger, " TOTAL globally frequent dfscodes = " << globally_frequent_dfscodes.size());

  external_mappings.clear();
  supports_locally_frequent.clear();
  supports_locally_not_frequent.clear();
  internal_supports.clear();
  external_supports.clear();

}

void global_support_handler::compute_globally_frequent_patterns_with_mappings(types::Graph &graph){

  timeval start_time, stop_time;
  map<string, double> profile_timings;
  string profiles[] = {"communicate locally freq patterns", "sort locally freq patterns", "identify external vmaps", "prepare arrays and communicate vmaps", "allreduce and compute support", "others"};

  for(int i = 0; i < 6; i++)
    profile_timings[profiles[i]] = 0.0;

  gettimeofday(&start_time, 0);

  //gather all locally frequent dfs codes
  std::vector<types::DFSCode> all_locally_frequent;

  get_all_locally_frequent_patterns(all_locally_frequent);

  gettimeofday(&stop_time, 0);
  profile_timings[profiles[0]] += utils::get_time_diff(start_time, stop_time);

  gettimeofday(&start_time, 0);

  DEBUG(*logger, "The locally frequent ONES = " << all_locally_frequent.size());

  //all ranks have the same order of the frequent patterns
  std::sort(all_locally_frequent.begin(), all_locally_frequent.end());

  gettimeofday(&stop_time, 0);
  profile_timings[profiles[1]] += utils::get_time_diff(start_time, stop_time);

  gettimeofday(&start_time, 0);

  //for 4 ranks and 2 partitions , groups (0,1) (2,3)
  MPI_Comm rank_group;
  int group_id = rank / num_partitions;
  MPI_Comm_split(MPI_COMM_WORLD, group_id, rank, &rank_group);

  int num_groups = numtasks / num_partitions;
  int patterns_per_group =  (int) ceil( all_locally_frequent.size() * 1.0 / num_groups );
  int start_index = group_id * patterns_per_group;
  int end_index = start_index + patterns_per_group - 1;
  if(end_index > all_locally_frequent.size() - 1)
    end_index = all_locally_frequent.size() - 1;

  DEBUG(*logger, "rank " << rank << " group id " << group_id << " total frequent " << all_locally_frequent.size() << " start_index " << start_index << " end_index " << end_index);

  const int MAX_PATTERNS = 200;

  //marks whether start_index + i pattern is a false positive or not, 1 for yes, 0 for no , also support value
  int *false_positive_marker = 0;
  int false_positive_marker_size = 2;

  gettimeofday(&stop_time, 0);
  profile_timings[profiles[5]] += utils::get_time_diff(start_time, stop_time);

  if(start_index <= end_index) {

    gettimeofday(&start_time, 0);

    false_positive_marker_size = 2 + (end_index - start_index + 1);
    false_positive_marker = new int[false_positive_marker_size];
    false_positive_marker[0] = start_index;
    false_positive_marker[1] = end_index;

    //handle patterns in parts
    int k = start_index;

    gettimeofday(&stop_time, 0);
    profile_timings[profiles[5]] += utils::get_time_diff(start_time, stop_time);

    while(k <= end_index) {

      gettimeofday(&start_time, 0);

      std::map<string, std::vector<std::vector<std::set<int> > > > dfscode_vertex_mappings;
      int num_patterns = 0, total_vertex_map_size = 0;

      if ( k + MAX_PATTERNS > end_index) {
        num_patterns = end_index - k + 1;
      }else{
        num_patterns = MAX_PATTERNS;
      }

      int *vertex_map_sizes = new int[num_patterns];

      //#pragma omp parallel for//reduction(+:total_vertex_map_size)
      //for(int i = 0; i < num_patterns; i++){
      numthreads = 1;
      #pragma omp parallel num_threads(numthreads)
      {
        int thread_id = omp_get_thread_num();
        int patterns_per_thread = (int) ceil(num_patterns * 1.0 / numthreads);
        int start_index2 = thread_id * patterns_per_thread;
        int end_index2 = start_index2 + patterns_per_thread - 1;
        if(end_index2 > num_patterns - 1)
          end_index2 = num_patterns - 1;

        for(int i = start_index2; i <= end_index2; i++) {

          types::DFSCode dfscode;
          //try{
          dfscode = all_locally_frequent.at(k + i);
          //}catch(std::out_of_range o){
          //	 std::cout<<o.what()<<std::endl;
          //}
          int vmap_size = 1;
          for(int m = 0; m <dfscode.size(); m++) {
            if(dfscode[m].is_forward())
              vmap_size++;
          }

          string dfs_str = dfscode.to_string();

          omp_set_lock(&lock_gs);
          std::vector<std::set<int> > vertex_maps(vmap_size);
          if(vertex_mappings.count(dfscode.to_string()) > 0)
            vertex_maps = vertex_mappings[dfscode.to_string()];
          omp_unset_lock(&lock_gs);

          //get_vertex_maps(dfscode, vertex_maps);

          std::vector<std::vector<std::set<int> > > global_ids_by_partitions(num_partitions, std::vector<std::set<int> >(vertex_maps.size() ));
          omp_set_lock(&lock_gs);
          dfscode_vertex_mappings[dfs_str] = global_ids_by_partitions;
          omp_unset_lock(&lock_gs);
          vertex_map_sizes[i] = vertex_maps.size();
          //total_vertex_map_size += vertex_maps.size();
          try{
            for(int j = 0; j < vertex_maps.size(); j++) {
              for( std::set<int>::iterator it = vertex_maps.at(j).begin(); it != vertex_maps.at(j).end(); it++) {
                int local_vid = *it;
                omp_set_lock(&lock_gs);
                //dfscode_vertex_mappings[dfs_str][graph[local_vid].orig_part_id][j].insert(graph[local_vid].global_vid);
                //dfscode_vertex_mappings.at(dfs_str).at(graph[local_vid].orig_part_id).at(j).insert(graph[local_vid].global_vid);
                dfscode_vertex_mappings[dfs_str].at(graph[local_vid].orig_part_id).at(j).insert(graph[local_vid].global_vid);
                DEBUG(*logger, " rank " << rank << " local vid " << local_vid << " owner rank = " << graph[local_vid].orig_part_id);
                omp_unset_lock(&lock_gs);

              }
            }
          } catch(std::out_of_range o) {
            std::cout << o.what() << std::endl;
          }
        }                         //for
      }                   // omp parallel
                          //}

      DEBUG(*logger, "Processing number of patterns " << num_patterns);

      for(int j = 0; j < num_patterns; j++)
        total_vertex_map_size += vertex_map_sizes[j];
      delete [] vertex_map_sizes;

      gettimeofday(&stop_time, 0);
      profile_timings[profiles[2]] += utils::get_time_diff(start_time, stop_time);

      gettimeofday(&start_time, 0);

      int *send_lengths = new int[num_partitions];           //number of neighbor requests per rank
      int *recv_lengths = new int[num_partitions];
      int total_send_length = 0, total_recv_length = 0;

      //#pragma omp parallel for
      for(int j = 0; j < num_partitions; j++) {
        send_lengths[j] = 0;
        if(j != rank % num_partitions) {
          for(int i = 0; i < num_patterns; i++) {
            types::DFSCode dfscode = all_locally_frequent[k + i];
            string dfs_str = dfscode.to_string();

            send_lengths[j] += dfscode_vertex_mappings[dfs_str][j].size();                                     // number of mappings for pattern vertices
            for(int l = 0; l < dfscode_vertex_mappings[dfs_str][j].size(); l++)
              send_lengths[j] += dfscode_vertex_mappings[dfs_str][j][l].size();
          }
          total_send_length += send_lengths[j];
        }
      }

      MPI_Alltoall(send_lengths, 1, MPI_INT, recv_lengths, 1, MPI_INT, rank_group);
      DEBUG(*logger," send lengths " << utils::print_array(send_lengths, num_partitions) << " recv lengths " << utils::print_array(recv_lengths, num_partitions));

      int *cum_recv_lengths = new int[num_partitions];
      //#pragma omp parallel for
      for(int j = 0; j < num_partitions; j++) {
        total_recv_length += recv_lengths[j];
      }

      //Now allocate send_buf and recv_buf;
      int *send_buf = new int[total_send_length];
      int *recv_buf = new int[total_recv_length];

      int offset = 0;
      for(int j = 0; j < num_partitions; j++) {
        if(j != rank % num_partitions) {
          for(int i = 0; i < num_patterns; i++) {
            types::DFSCode dfscode = all_locally_frequent[k + i];
            string dfs_str = dfscode.to_string();

            for(int l = 0; l < dfscode_vertex_mappings[dfs_str][j].size(); l++) {
              send_buf[offset++] = dfscode_vertex_mappings[dfs_str][j][l].size();
              for(std::set<int>::iterator it  = dfscode_vertex_mappings[dfs_str][j][l].begin(); it!= dfscode_vertex_mappings[dfs_str][j][l].end(); it++)
                send_buf[offset++] =  *it;
            }
          }
        }
      }

      //Perform MPI_Alltoallv
      int *sdispls = new int[num_partitions];
      int *rdispls = new int[num_partitions];
      sdispls[0] = 0;
      rdispls[0] = 0;

      for(int j = 1; j< num_partitions; j++) {
        sdispls[j] = sdispls[j - 1] + send_lengths[j - 1];
        rdispls[j] = rdispls[j - 1] + recv_lengths[j - 1];
      }

      MPI_Alltoallv(send_buf, send_lengths, sdispls, MPI_INT, recv_buf, recv_lengths, rdispls, MPI_INT, rank_group);
      delete [] rdispls;
      delete [] sdispls;

      DEBUG(*logger, " send_buf " << utils::print_array(send_buf, total_send_length) << " recv_buf " << utils::print_array(recv_buf, total_recv_length));

      gettimeofday(&stop_time, 0);
      profile_timings[profiles[3]] += utils::get_time_diff(start_time, stop_time);

      gettimeofday(&start_time, 0);

      //Now process received data, compute actual local support
      offset = 0;
      for(int j = 0; j < num_partitions; j++) {
        if(j != rank % num_partitions) {
          for(int i = 0; i < num_patterns; i++) {
            types::DFSCode dfscode = all_locally_frequent[k + i];
            string dfs_str = dfscode.to_string();

            for(int l = 0; l < dfscode_vertex_mappings[dfs_str][rank % num_partitions].size(); l++) {
              int len = recv_buf[offset++];
              //omp_set_lock(&lock_gs);
              for(int m = 0; m < len; m++) {
                dfscode_vertex_mappings[dfs_str][rank % num_partitions][l].insert(recv_buf[offset++]);
              }
              //omp_unset_lock(&lock_gs);
            }
          }
        }
      }

      //Allreduce
      int *sendcounts = new int[total_vertex_map_size];
      int *recvcounts = new int[total_vertex_map_size];

      offset = 0;
      for(int i = 0; i < num_patterns; i++) {
        types::DFSCode dfscode = all_locally_frequent[k + i];
        string dfs_str = dfscode.to_string();
        for(int l = 0; l < dfscode_vertex_mappings[dfs_str][rank % num_partitions].size(); l++) {                        //vertex map size
          sendcounts[offset++] = dfscode_vertex_mappings[dfs_str][rank % num_partitions][l].size();
        }
      }

      MPI_Allreduce(sendcounts, recvcounts, total_vertex_map_size, MPI_INT, MPI_SUM, rank_group);

      //compute support
      offset = 0;

      for(int i = 0; i < num_patterns; i++) {
        types::DFSCode dfscode = all_locally_frequent[k + i];
        string dfs_str = dfscode.to_string();
        int min = 0x7FFFFFFF;
        for(int l = 0; l < dfscode_vertex_mappings[dfs_str][rank % num_partitions].size(); l++) {                        //vertex map size
          if(min > recvcounts[offset])
            min = recvcounts[offset];
          offset++;
        }

        DEBUG(*logger, "DFS code " << dfscode.to_string() << " vertex mappings " << utils::print_array(recvcounts + (offset - 3),3));

        if(min < minsupp) {
          //if(rank == 0)
          //num_false_positives++;
          DEBUG(*logger, "DFS code" << dfscode.to_string() << " has support " << min << " is not globally frequent");
          //Remove if present in local list
          //remove_globally_frequent_dfscode(dfscode);
          false_positive_marker[2 + k - start_index + i] = 1;
        }else{

          if(supports_locally_frequent.count(dfscode.to_string()) > 0 ||
             supports_locally_not_frequent.count(dfscode.to_string()) > 0) {
            globally_frequent_dfscodes.push_back(dfscode);
          }

          //all_globally_estimated_dfscodes.push_back(dfscode);
          //update current support
          //set_dfscode_support(dfscode, min);
          false_positive_marker[2 + k - start_index + i] = 0;
          if(rank == 0) {
            //report(dfscode, min);
            omp_set_lock(&lock_gs);
            frequent_patterns++;
            omp_unset_lock(&lock_gs);
            DEBUG(*logger, "Globally frequent = " << dfscode.to_string() << " : support = " << min);
          }
        }
      }

      gettimeofday(&stop_time, 0);
      profile_timings[profiles[4]] += utils::get_time_diff(start_time, stop_time);

      gettimeofday(&start_time, 0);

      delete [] send_lengths;
      delete [] recv_lengths;
      delete [] send_buf;
      delete [] recv_buf;
      delete [] sendcounts;
      delete [] recvcounts;
      //delete [] send_buf2;
      //delete [] recv_buf2;
      delete [] cum_recv_lengths;

      k += MAX_PATTERNS;


      gettimeofday(&stop_time, 0);
      profile_timings[profiles[5]] += utils::get_time_diff(start_time, stop_time);

    }

  }else{
    false_positive_marker = new int[2];
    false_positive_marker[0] = start_index;
    false_positive_marker[1] = end_index;
    false_positive_marker_size = 2;
    //do nothing
  }

  if(rank == 0)
    INFO(*logger, "GLOBAL SUPPORT: level " << dfscode_size << " total locally frequent patterns " << all_locally_frequent.size());


  for(int i = 0; i < 6; i++) {
    double sum = 0.0;
    MPI_Reduce(&profile_timings[profiles[i]], &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //cout<<"rank "<<rank<<" profile = "<<profiles[i]<<" time = "<<profile_timings[profiles[i]]<<endl;
    if(rank == 0)
      INFO(*logger, "GLOBAL SUPPORT: level " << dfscode_size << " AVG time for " << profiles[i] << " = " << (sum / numtasks));
  }

  //now send info to the mirror/peer ranks
  //Let the mirror processor know of the frequent patterns
  //For example for 2 partitions and 4 processors (0,2) and (1,3) should share
  //independently computed info among them

  if(numtasks > num_partitions) {
    // communicate to all ranks handling the same partition
    MPI_Comm rank_mirror;
    int group_id2 = rank % num_partitions;
    MPI_Comm_split(MPI_COMM_WORLD, group_id2, rank, &rank_mirror);

    int num_mirrors = numtasks / num_partitions;
    int *lens = new int[num_mirrors];

    MPI_Allgather(&false_positive_marker_size, 1, MPI_INT, lens, 1, MPI_INT, rank_mirror);

    int total_recv_length = 0;
    for(int i = 0; i < num_mirrors; i++)
      total_recv_length += lens[i];

    int *recv_buf2 = new int[total_recv_length];

    DEBUG(*logger, " Lengths received: " << utils::print_array(lens, num_mirrors));

    int *displs = new int[num_mirrors];
    displs[0] = 0;
    for(int i = 1; i< num_mirrors; i++)
      displs[i] = displs[i - 1] + lens[i - 1];

    MPI_Allgatherv(false_positive_marker, false_positive_marker_size, MPI_INT, recv_buf2, lens, displs, MPI_INT, rank_mirror);
    delete [] displs;

    //INFO(*logger, " Sent "<<utils::print_array(false_positive_marker, false_positive_marker_size) );
    //INFO(*logger, " Received "<<utils::print_array(recv_buf2, total_recv_length));

    //now process the false positives
    int offset = 0;
    for(int i = 0; i < num_mirrors; i++) {
      if( rank / num_partitions == i ) {
        //already processed, advance offset
        offset += lens[i];
      }else{

        int si = recv_buf2[offset++];
        int ei = recv_buf2[offset++];

        if(si <= ei) {
          for(int j = si; j <= ei; j++) {
            types::DFSCode dfscode = all_locally_frequent[j];
            int is_false_pos = recv_buf2[offset++];
            //int sup = recv_buf2[offset++];
            DEBUG(*logger, "DFS Code " << all_locally_frequent[j].to_string() << " is false pos" << is_false_pos << " sup " << sup);
            if(is_false_pos) {
              //if(rank == 0)
              //num_false_positives++;
              //Remove if present in local list
              remove_globally_frequent_dfscode(dfscode);
            }else{
              //update current support
              //gs_handler.set_dfscode_support(dfscode, sup);
              DEBUG(*logger, "DFS code" << dfscode.to_string() << " has support " << min << " is not globally frequent");
            }
          }
        }

      }
    }

    delete [] lens;
    delete [] recv_buf2;
  }

  delete [] false_positive_marker;

}

void global_support_handler::get_all_locally_frequent_patterns_with_internal_supports(std::vector<types::DFSCode> &all_locally_frequent){

  std::set<string> locally_frequent_set;

  int my_num_patterns = supports_locally_frequent.size();
  int *all_num_patterns = new int[numtasks];
  all_gather_msg(&my_num_patterns, all_num_patterns, 1);

  int *recv_lengths = new int[numtasks];
  int total_recv_len = 0;
  stringstream ss;
  ss << " number of globally frequent patterns = ";
  for(int i = 0; i<numtasks; i++) {
    total_recv_len += all_num_patterns[i] * dfscode_size * 5;
    recv_lengths[i] = all_num_patterns[i] * dfscode_size * 5;
    ss << all_num_patterns[i] << " ";
  }

  DEBUG(*logger, ss.str());

  //int *send_buf = new int[my_num_patterns * ( dfscode_size * 5 + 1 + (dfscode_size + 1) ) ];
  int *send_buf = new int[my_num_patterns * dfscode_size * 5];
  int *recv_buf = new int[total_recv_len];

  int offset = 0;
  //cout<<"rank " << rank <<" sending DFS codes ";
  for(std::map<string, std::vector<int> >::iterator it = supports_locally_frequent.begin(); it != supports_locally_frequent.end(); it++) {
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(it->first);
    for(int i = 0; i < DFS_CODE.size(); i++) {
      DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
      offset += 5;
    }

    locally_frequent_set.insert(it->first);

  }

  all_gatherv_msg(send_buf, recv_buf, my_num_patterns * dfscode_size * 5, recv_lengths);

  ss.clear();
  ss << "rank " << rank << " Received locally frequent dfscodes ";
  offset = 0;
  for(int k = 0; k<numtasks; k++) {
    for(int i = 0; i<recv_lengths[k]; i++)
      ss << recv_buf[offset + i] <<  " ";
    ss << " | ";
    offset += recv_lengths[k];
  }
  DEBUG(*logger, ss.str());

  offset = 0;
  for(int i = 0; i<numtasks; i++) {
    if(i != rank) {

      int *buffer = recv_buf + offset;
      int recv_pos = 0;
      for(int k = 0; k < all_num_patterns[i]; k++) {
        //get DFS code
        types::DFSCode DFS_CODE;
        for(int j = 0; j < dfscode_size; j++) {
          types::DFS dfs;
          dfs.deserialize((char*) (buffer + recv_pos), 5 * sizeof(int));
          DFS_CODE.push_back(dfs);
          recv_pos += 5;
        }

        string dfscode = DFS_CODE.to_string();
        DEBUG(*logger, "Locally frequent dfscode received " << DFS_CODE.to_string());
        locally_frequent_set.insert(dfscode);

      }
    }
    offset += recv_lengths[i];
  }

  int send_support_length = 0, recv_support_length = 0;
  for(std::set<string>::iterator it = locally_frequent_set.begin(); it != locally_frequent_set.end(); it++) {
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(*it);
    all_locally_frequent.push_back(DFS_CODE);
    send_support_length++;             // support size 0 if does not exist in this partition
    if(internal_supports.count(*it) > 0)
      send_support_length += internal_supports[*it].size();

  }

  MPI_Allreduce(&send_support_length, &recv_support_length, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  all_gather_msg(&send_support_length, recv_lengths, 1);

  delete [] send_buf;
  delete [] recv_buf;
  send_buf = new int[send_support_length];
  recv_buf = new int[recv_support_length];

  offset = 0;
  for(std::set<string>::iterator it = locally_frequent_set.begin(); it != locally_frequent_set.end(); it++) {
    types::DFSCode DFS_CODE = types::DFSCode::read_from_str(*it);
    DEBUG(*logger, "DFS CODE = " << *it << " internal support " << utils::print_vector(global_supports_internal[*it]));
    if(internal_supports.count(*it) > 0 ) {
      send_buf[offset++] = internal_supports[*it].size();
      for(int j = 0; j < internal_supports[*it].size(); j++)
        send_buf[offset++] = internal_supports[*it][j];
      global_supports[*it] = internal_supports[*it];
      global_supports_internal[*it] = internal_supports[*it];
    } else{
      send_buf[offset++] = 0;
    }

  }

  all_gatherv_msg(send_buf, recv_buf, send_support_length, recv_lengths);

  DEBUG(*logger, "RECV LENGTHS = " << utils::print_array(recv_lengths,numtasks) << " recv_buf = " << utils::print_array(recv_buf,recv_support_length));

  int start_offset = 0;
  for(int i = 0; i < numtasks; i++) {
    offset = start_offset;
    if (i != rank) {
      for(std::set<string>::iterator it = locally_frequent_set.begin(); it != locally_frequent_set.end(); it++) {
        int supp_len = recv_buf[offset++];
        if(supp_len > 0) {
          if(global_supports.count(*it) == 0) {
            std::vector<int> t;
            global_supports[*it] = t;
            global_supports_internal[*it] = t;
            for(int j = 0; j < supp_len; j++) {
              global_supports[*it].push_back(0);
              global_supports_internal[*it].push_back(0);
            }
          }

          for(int j = 0; j <supp_len; j++) {
            global_supports[*it][j] += recv_buf[offset];
            global_supports_internal[*it][j] += recv_buf[offset++];
          }
        }
      }
    }
    start_offset += recv_lengths[i];
  }

  DEBUG(*logger, "Total locally frequent = " << all_locally_frequent.size());
  delete [] all_num_patterns;
  delete [] send_buf;
  delete [] recv_buf;
  delete [] recv_lengths;

}


void global_support_handler::get_remaining_frequent_patterns(std::vector<types::DFSCode> &remaining_estimated){

  //std::set<string> frequent_set;

  int my_num_patterns = 0;
  for(std::vector<types::DFSCode>::iterator it = globally_frequent_dfscodes_hashed.begin(); it != globally_frequent_dfscodes_hashed.end(); it++) {
    if(marked_actual_globally_frequent_patterns[it->to_string()] == 0) {            //not actually frequent, not a candidate for false pos
      my_num_patterns++;
    }
  }

  int *all_num_patterns = new int[numtasks];
  all_gather_msg(&my_num_patterns, all_num_patterns, 1);

  int *recv_lengths = new int[numtasks];
  int total_recv_len = 0;
  stringstream ss;
  ss << " number of globally frequent patterns = ";
  for(int i = 0; i<numtasks; i++) {
    total_recv_len += all_num_patterns[i] * dfscode_size * 5;
    recv_lengths[i] = all_num_patterns[i] * dfscode_size * 5;
    ss << all_num_patterns[i] << " ";
  }

  DEBUG(*logger, ss.str());

  int *send_buf = new int[my_num_patterns * dfscode_size * 5];
  int *recv_buf = new int[total_recv_len];

  int offset = 0;
  //cout<<"rank " << rank <<" sending DFS codes ";
  for(std::vector<types::DFSCode>::iterator it = globally_frequent_dfscodes_hashed.begin(); it != globally_frequent_dfscodes_hashed.end(); it++) {
    if(marked_actual_globally_frequent_patterns[it->to_string()] == 1) {            //already actually frequent, not a candidate for false pos
      continue;
    }

    types::DFSCode DFS_CODE = *it;             //types::DFSCode::read_from_str(it->first);
    for(int i = 0; i < DFS_CODE.size(); i++) {
      DFS_CODE[i].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
      offset += 5;
    }

    //frequent_set.insert(it->to_string());
    //remaining_estimated.push_back(*it);
  }
  //cout<<endl;

  //all_gatherv_msg(send_buf, recv_buf, my_num_patterns * ( dfscode_size * 5 + 1 + (dfscode_size + 1) ), recv_lengths);
  all_gatherv_msg(send_buf, recv_buf, my_num_patterns * dfscode_size * 5, recv_lengths);

  ss.clear();
  ss << "rank " << rank << " Received locally frequent dfscodes ";
  offset = 0;
  for(int k = 0; k<numtasks; k++) {
    for(int i = 0; i<recv_lengths[k]; i++)
      ss << recv_buf[offset + i] <<  " ";
    ss << " | ";
    offset += recv_lengths[k];
  }
  DEBUG(*logger, ss.str());

  offset = 0;
  for(int i = 0; i<numtasks; i++) {
    //if(i != rank){
    //process_globally_frequent_patterns(recv_buf + offset, recv_lengths[i], all_num_patterns[i], i);
    int *buffer = recv_buf + offset;
    //int len = recv_lengths[i];
    int recv_pos = 0;
    //while(recv_pos < len){
    for(int k = 0; k < all_num_patterns[i]; k++) {
      //get DFS code
      types::DFSCode DFS_CODE;
      for(int j = 0; j < dfscode_size; j++) {
        types::DFS dfs;
        dfs.deserialize((char*) (buffer + recv_pos), 5 * sizeof(int));
        DFS_CODE.push_back(dfs);
        recv_pos += 5;
      }

      string dfscode = DFS_CODE.to_string();
      DEBUG(*logger, "remaining frequent dfscode received " << DFS_CODE.to_string());
      //frequent_set.insert(dfscode);
      remaining_estimated.push_back(DFS_CODE);
    }
    //}
    offset += recv_lengths[i];
  }

  DEBUG(*logger, "Total remaining frequent = " << remaining_estimated.size());
  delete [] all_num_patterns;
  delete [] send_buf;
  delete [] recv_buf;
  delete [] recv_lengths;
}



int global_support_handler::process_non_frequent_local_support_requests(int *buffer, int len, int num_dfscodes, int dest, int *data, int &datasize){

  int found_dfscodes = 0, received_dfscodes = 0;
  int recv_pos = 0;
  int send_pos = 0;

  while(recv_pos < len) {
    //get DFS code
    types::DFSCode DFS_CODE;
    for(int i = 0; i < dfscode_size; i++) {
      types::DFS dfs;
      dfs.deserialize((char*) (buffer + recv_pos), 5 * sizeof(int));
      DFS_CODE.push_back(dfs);
      recv_pos += 5;
    }

    DEBUG(*logger, "DFS code requested " << DFS_CODE.to_string());
    if(supports_locally_not_frequent.count(DFS_CODE.to_string()) > 0) {
      int support_size = supports_locally_not_frequent[DFS_CODE.to_string()].size();
      data[send_pos++] = support_size;
      for(int i = 0; i < support_size; i++)
        data[send_pos++] = supports_locally_not_frequent[DFS_CODE.to_string()][i];
      found_dfscodes++;
    }
    else{
      data[send_pos++] = -1;                   //not found, invalid entry for the dfs code
    }
    received_dfscodes++;
  }

  datasize += send_pos;

  DEBUG(*logger, "Rank " << rank << " processed dfscode requests " << received_dfscodes << " from " << dest);
  return found_dfscodes;
}

int global_support_handler::process_non_frequent_local_support_requests_with_int_ext_map_cardinalities(int *buffer, int len, int num_dfscodes, int dest, int *data, int &datasize){

  int found_dfscodes = 0, received_dfscodes = 0;
  int recv_pos = 0;
  int send_pos = 0;

  while(recv_pos < len) {
    //get DFS code
    types::DFSCode DFS_CODE;
    for(int i = 0; i < dfscode_size; i++) {
      types::DFS dfs;
      dfs.deserialize((char*) (buffer + recv_pos), 5 * sizeof(int));
      DFS_CODE.push_back(dfs);
      recv_pos += 5;
    }

    DEBUG(*logger, "DFS code requested " << DFS_CODE.to_string());
    if(supports_locally_not_frequent.count(DFS_CODE.to_string()) > 0) {
      int support_size = supports_locally_not_frequent[DFS_CODE.to_string()].size();
      data[send_pos++] = support_size;
      for(int i = 0; i < support_size; i++) {
        data[send_pos++] = internal_supports[DFS_CODE.to_string()][i];
        data[send_pos++] = external_supports[DFS_CODE.to_string()][i];
      }
      found_dfscodes++;
    }
    else{
      data[send_pos++] = -1;                   //not found, invalid entry for the dfs code
    }
    received_dfscodes++;
  }

  datasize += send_pos;

  DEBUG(*logger, "Rank " << rank << " processed dfscode requests " << received_dfscodes << " from " << dest);
  return found_dfscodes;
}


int global_support_handler::process_non_frequent_local_support_requests(int *buffer, int len, int num_dfscodes, int dest){

  int found_dfscodes = 0;
  int recv_pos = 0;
  int send_pos = 0;
  //allocate max possible size for data; max vertices = dfscode_size + 1
  int send_size = 1 /*MSG_SUPPORT_DATA*/ + num_dfscodes * ( 1 /*support size*/ + dfscode_size * 5  + (dfscode_size + 1) /*max num vertices*/ );
  //int *data = new int[send_size];
  int *data;
  NEW_INT_ARRAY(data, send_size, *logger);
  data[send_pos++] = MSG_SUPPORT_DATA;
  while(recv_pos < len) {
    //get DFS code
    types::DFSCode DFS_CODE;
    for(int i = 0; i < dfscode_size; i++) {
      types::DFS dfs;
      dfs.deserialize((char*) (buffer + recv_pos), 5 * sizeof(int));
      DFS_CODE.push_back(dfs);
      recv_pos += 5;
    }

    DEBUG(*logger, "DFS code requested " << DFS_CODE.to_string());
    if(supports_locally_not_frequent.count(DFS_CODE.to_string()) > 0) {
      int support_size = supports_locally_not_frequent[DFS_CODE.to_string()].size();
      data[send_pos++] = support_size;
      for(int i = 0; i < DFS_CODE.size(); i++) {
        DFS_CODE[i].serialize( (char*) (data + send_pos), sizeof(int) * 5 );
        send_pos += 5;
      }

      for(int i = 0; i < support_size; i++)
        data[send_pos++] = supports_locally_not_frequent[DFS_CODE.to_string()][i];

      found_dfscodes++;

      //does not matter anyway, it is not requested since locally not frequent
      //supports_locally_not_frequent.erase(DFS_CODE.to_string());
    }
  }

  if(found_dfscodes > 0) {
    //DEBUG(*logger, "msg ID " << dfs_id << " DFS code support found " << DFS_CODE.to_string());
    send_data_msg(data, send_pos, dest);
  }else{
    //DEBUG(*logger, "msg ID " << dfs_id << " DFS code support not found " << DFS_CODE.to_string());
  }

  //delete []data;
  DELETE_INT_ARRAY(data, *logger);
  return found_dfscodes;
}

void global_support_handler::process_globally_frequent_patterns(int *buffer, int len, int num_dfscodes, int dest){

  int recv_pos = 0;
  while(recv_pos < len) {
    //get DFS code
    types::DFSCode DFS_CODE;
    for(int i = 0; i < dfscode_size; i++) {
      types::DFS dfs;
      dfs.deserialize((char*) (buffer + recv_pos), 5 * sizeof(int));
      DFS_CODE.push_back(dfs);
      recv_pos += 5;
    }

    DEBUG(*logger, "Globally frequent dfscode received " << DFS_CODE.to_string());
    if(supports_locally_frequent.count(DFS_CODE.to_string()) > 0 ||
       supports_locally_not_frequent.count(DFS_CODE.to_string()) > 0) {
      globally_frequent_dfscodes.push_back(DFS_CODE);
    }

    all_globally_estimated_dfscodes.push_back(DFS_CODE);
  }

}

void global_support_handler::process_globally_frequent_patterns_with_frequent_flags(int *buffer, int len, int num_dfscodes, int dest){

  int recv_pos = 0;
  while(recv_pos < len) {
    //get DFS code
    types::DFSCode DFS_CODE;
    for(int i = 0; i < dfscode_size; i++) {
      types::DFS dfs;
      dfs.deserialize((char*) (buffer + recv_pos), 5 * sizeof(int));
      DFS_CODE.push_back(dfs);
      recv_pos += 5;
    }

    int freq_flag = buffer[recv_pos++];

    DEBUG(*logger, "Globally frequent dfscode received " << DFS_CODE.to_string());
    if(supports_locally_frequent.count(DFS_CODE.to_string()) > 0 ||
       supports_locally_not_frequent.count(DFS_CODE.to_string()) > 0) {
      globally_frequent_dfscodes.push_back(DFS_CODE);
    }

    if(freq_flag == 0)
      all_globally_estimated_dfscodes.push_back(DFS_CODE);
  }

}

void global_support_handler::process_globally_frequent_patterns2(int *buffer, int len, int num_dfscodes){

  int recv_pos = 0;
  while(recv_pos < len) {
    //get DFS code
    types::DFSCode DFS_CODE;
    for(int i = 0; i < dfscode_size; i++) {
      types::DFS dfs;
      dfs.deserialize((char*) (buffer + recv_pos), 5 * sizeof(int));
      DFS_CODE.push_back(dfs);
      recv_pos += 5;
    }
    globally_frequent_dfscodes.push_back(DFS_CODE);
  }
}



void global_support_handler::process_received_data(int *&recv_buf, int len, int source){

  int pos = 1;       // skip MSG_SUPPORT_DATA

  int dfscodes_received = 0;

  //if(rank == 1)
  DEBUG(*logger, "Rank " << rank << " received data " << utils::print_array(recv_buf, len));

  while(pos < len) {
    int support_size = recv_buf[pos++];

    types::DFSCode DFS_CODE;

    for(int i = 0; i < dfscode_size; i++) {
      types::DFS dfs;
      dfs.deserialize((char*) (recv_buf + pos), 5 * sizeof(int));
      DFS_CODE.push_back(dfs);
      pos += 5;
    }

    string dfscode_str = DFS_CODE.to_string();

    DEBUG(*logger, "rank " << rank << " received hashed dfscode = " << dfscode_str);

    if(global_supports.count(dfscode_str) == 0) {
      std::vector<int> t;
      global_supports[dfscode_str] = t;

      for(int i = 0; i <support_size; i++)
        global_supports[dfscode_str].push_back(recv_buf[pos++]);

      //if could be in locally not frequent list
      if(supports_locally_not_frequent.count(dfscode_str) > 0) {
        for(int i = 0; i <support_size; i++)
          global_supports[dfscode_str][i] += supports_locally_not_frequent[dfscode_str][i];
      }

    } else {

      for(int i = 0; i <support_size; i++)
        global_supports[dfscode_str][i] += recv_buf[pos++];
    }

    dfscodes_received++;
    patterns_received++;
  }

  DEBUG(*logger, "rank " << rank << " Total dfscodes received from " << source << " = " << dfscodes_received << " " << patterns_received);
}


void global_support_handler::process_received_data2(int *recv_buf, int len, int source){

  int pos = 0;
  int dfscodes_received = 0;

  //if(rank == 1)
  DEBUG(*logger, "Rank " << rank << " received data from " << source << " = " << utils::print_array(recv_buf, len));

  for(std::map<string, std::vector<int> >::iterator it = global_supports.begin(); it != global_supports.end(); ++it) {

    dfscodes_received++;

    if(recv_buf[pos] == -1) {
      pos++;
      continue;
    }

    int support_size = recv_buf[pos++];
    DEBUG(*logger, " global supports " << utils::print_vector(it->second));
    for(int i = 0; i <support_size; i++)
      (it->second)[i] += recv_buf[pos++];
  }

  DEBUG(*logger, "rank " << rank << " *** Total dfscodes received from " << source << " = " << dfscodes_received );
}

void global_support_handler::process_received_data2_with_int_ext_map_cardinalities(int *recv_buf, int len, int source){

  int pos = 0;
  int dfscodes_received = 0;

  //if(rank == 1)
  DEBUG(*logger, "Rank " << rank << " received data from " << source << " = " << utils::print_array(recv_buf, len));

  for(std::map<string, std::vector<int> >::iterator it = global_supports.begin(); it != global_supports.end(); ++it) {

    dfscodes_received++;

    if(recv_buf[pos] == -1) {
      pos++;
      continue;
    }

    int support_size = recv_buf[pos++];
    DEBUG(*logger, " global supports " << utils::print_vector(it->second));
    for(int i = 0; i <support_size; i++) {
      global_supports_internal[it->first][i] += recv_buf[pos++];
      global_supports_external[it->first][i] += recv_buf[pos++];
      (it->second)[i] += recv_buf[pos - 2] + recv_buf[pos - 1];
      if(global_supports_external[it->first][i] < recv_buf[pos - 1])
        global_supports_external[it->first][i] = recv_buf[pos - 1];
    }
  }

  DEBUG(*logger, "rank " << rank << " *** Total dfscodes received from " << source << " = " << dfscodes_received );
}


int global_support_handler::get_globally_frequent(){
  return frequent_patterns;
}

/*void global_support_handler::clear_dfscode_ids(){
        dfscode_sent_ids.clear();
        dfscode_received_ids.clear();
        received_local_support_from_ids.clear();
        globally_frequent_ids.clear();
        global_supports.clear();
   }*/

void global_support_handler::clear_vertex_mappings(){
  vertex_mappings.clear();
  external_mappings.clear();
  vertex_mappings_file_pos.clear();
  vertex_mappings_size_in_file.clear();
}

int global_support_handler::get_total_frequent_dfscodes(){
  int sendbuf = globally_frequent_dfscodes.size(), recvbuf = 0;

  MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  DEBUG(*logger, " sent " << sendbuf << " recv " << recvbuf);
  return recvbuf;

}

void global_support_handler::populate_globally_frequent_dfscodes(std::vector<types::DFSCode> &dfscodes_to_process){

  for(int i = 0; i<globally_frequent_dfscodes.size(); i++) {
    dfscodes_to_process.push_back(globally_frequent_dfscodes[i]);
  }
  globally_frequent_dfscodes.clear();
  all_globally_estimated_dfscodes.clear();
}

void global_support_handler::print_average_vertex_mappings_info(){

  std::vector<int> total_unique_vertex_mappings(dfscode_size + 1);       //allocated max size;
  std::vector<int> contributing_patterns_count(dfscode_size + 1);
  std::vector<int> average_unique_vertex_mappings(dfscode_size + 1);       //allocated max size;

  for(int i = 0; i <dfscode_size + 1; i++) {
    total_unique_vertex_mappings[i] = 0;
    contributing_patterns_count[i] = 0;
  }

  //for locally frequent
  for(std::map<string, vector<int> >::iterator it = supports_locally_frequent.begin(); it != supports_locally_frequent.end(); it++) {
    for(int i = 0; i< (it->second).size(); i++) {
      total_unique_vertex_mappings[i] += (it->second)[i];
      contributing_patterns_count[i]++;
    }
  }

  for(int i = 0; i<dfscode_size + 1; i++)
    average_unique_vertex_mappings[i] = total_unique_vertex_mappings[i] * 1.0 / contributing_patterns_count[i];
  INFO(*logger, "rank = " << rank << " level = " << dfscode_size << " locally frequent patterns vertex mappings: total = " << utils::print_vector(total_unique_vertex_mappings) << " : average = " << utils::print_vector(average_unique_vertex_mappings));

  //reset for locally not frequent
  for(int i = 0; i <dfscode_size + 1; i++) {
    total_unique_vertex_mappings[i] = 0;
    contributing_patterns_count[i] = 0;
  }

  for(std::map<string, vector<int> >::iterator it = supports_locally_not_frequent.begin(); it != supports_locally_not_frequent.end(); it++) {
    for(int i = 0; i< (it->second).size(); i++) {
      total_unique_vertex_mappings[i] += (it->second)[i];
      contributing_patterns_count[i]++;
    }
  }

  for(int i = 0; i<dfscode_size + 1; i++)
    average_unique_vertex_mappings[i] = total_unique_vertex_mappings[i] * 1.0 / contributing_patterns_count[i];

  INFO(*logger, "rank = " << rank << " level = " << dfscode_size << " locally NOT frequent patterns vertex mappings: total = " << utils::print_vector(total_unique_vertex_mappings) << " : average = " << utils::print_vector(average_unique_vertex_mappings));

}

void global_support_handler::print_globally_frequent_dfscodes(){
  if(out) {
    for(int i = 0; i < globally_frequent_dfscodes.size(); i++) {
      //out->output_graph(globally_frequent_dfscodes[i], final_dfscode_supports[globally_frequent_dfscodes[i].to_string()]);
      out->output_graph(globally_frequent_dfscodes[i]);
    }
  }
}

void global_support_handler::print_profile_info(){
  INFO(*logger, " hash sendrecv time = " << hash_sendrecv_time
       << " compute support time = " << compute_support_time
       << " unfinished pattern sendrecv time = " << unfinished_sendrecv_time
       <<  " global pattern send recv time = " << global_pattern_sendrecv_time);
}

void global_support_handler::bcast_msg(int *buffer, int length, int root) {
  //DEBUG(*logger, dfscodes_sent << " BROADCASTING buffer of length " << length << " root " << root);
  MPI_Bcast(buffer, length, MPI_INT, root, MPI_COMM_WORLD);
}

void global_support_handler::all_gather_msg(int *sendbuf, int *recvbuf, int length) {
  //DEBUG(*logger, dfscodes_sent << " ALL_GATHER buffer of length " << length);
  MPI_Allgather(sendbuf, length, MPI_INT, recvbuf, length, MPI_INT, MPI_COMM_WORLD);
}

void global_support_handler::all_gather_msg(int *sendbuf, int *recvbuf, int length, MPI_Comm &rank_group) {
  MPI_Allgather(sendbuf, length, MPI_INT, recvbuf, length, MPI_INT, rank_group);
}

void global_support_handler::all_gatherv_msg(int *sendbuf, int *recvbuf, int sendlen, int *recvlen) {
  //DEBUG(*logger, dfscodes_sent << " ALL_GATHERV");
  int *displs = new int[numtasks];
  displs[0] = 0;
  for(int i = 1; i< numtasks; i++)
    displs[i] = displs[i - 1] + recvlen[i - 1];
  MPI_Allgatherv(sendbuf, sendlen, MPI_INT, recvbuf, recvlen, displs, MPI_INT, MPI_COMM_WORLD);
  delete [] displs;
}

void global_support_handler::all_gatherv_msg(int *sendbuf, int *recvbuf, int sendlen, int *recvlen, MPI_Comm &rank_group) {
  //DEBUG(*logger, dfscodes_sent << " ALL_GATHERV");
  int *displs = new int[numtasks];
  displs[0] = 0;
  for(int i = 1; i< numtasks; i++)
    displs[i] = displs[i - 1] + recvlen[i - 1];
  MPI_Allgatherv(sendbuf, sendlen, MPI_INT, recvbuf, recvlen, displs, MPI_INT, rank_group);
  delete [] displs;
}

void global_support_handler::alltoall_msg(int *sendbuf, int *recvbuf, int length) {
  //DEBUG(*logger, dfscodes_sent << " XXX ALLTOALL buffer of length " << length << " root " << root);
  MPI_Alltoall(sendbuf, length, MPI_INT, recvbuf, length, MPI_INT, MPI_COMM_WORLD);
}


void global_support_handler::alltoallv_msg(int *sendbuf, int *recvbuf, int *sendlen, int *recvlen, int *sdispls, int *rdispls) {
  //DEBUG(*logger, dfscodes_sent << " XXX ALLTOALLV buffer of length " << length << " root " << root);
  MPI_Alltoallv(sendbuf, sendlen, sdispls, MPI_INT, recvbuf, recvlen, rdispls, MPI_INT, MPI_COMM_WORLD);
}


void global_support_handler::alltoallv_msg(int *sendbuf, int *recvbuf, int *sendlen, int *recvlen) {
  //DEBUG(*logger, dfscodes_sent << " XXX ALLTOALLV buffer of length " << length << " root " << root);
  int *sdispls = new int[numtasks];
  int *rdispls = new int[numtasks];
  sdispls[0] = 0;
  rdispls[0] = 0;
  for(int i = 1; i< numtasks; i++) {
    sdispls[i] = sdispls[i - 1] + sendlen[i - 1];
    rdispls[i] = rdispls[i - 1] + recvlen[i - 1];
  }

  MPI_Alltoallv(sendbuf, sendlen, sdispls, MPI_INT, recvbuf, recvlen, rdispls, MPI_INT, MPI_COMM_WORLD);
  delete [] rdispls;
  delete [] sdispls;

}

void global_support_handler::send_data_msg(int *buffer, int length, int dest_proc) {
  //DEBUG(*logger, " XXX SENDING buffer of length " << length << " of type MSG_GLOBAL_SUPPORT to " << dest_proc);
  MPI_Isend(buffer, length, MPI_INT, dest_proc, MSG_GLOBAL_SUPPORT, MPI_COMM_WORLD, &request);
}

void global_support_handler::recv_data_msg(int *buffer, int length, int src_proc, MPI_Status &stat) {
  //DEBUG(*logger, " XXX RECEIVING buffer of length " << length << " of type MSG_GLOBAL_SUPPORT from " << src_proc);
  MPI_Recv(buffer, length, MPI_INT, src_proc, MSG_GLOBAL_SUPPORT, MPI_COMM_WORLD, &stat);
}

} //namespace alg
