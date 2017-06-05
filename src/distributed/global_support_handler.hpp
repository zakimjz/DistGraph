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

#ifndef __GLOBAL_SUPPORT_HANDLER_HPP__
#define __GLOBAL_SUPPORT_HANDLER_HPP__

#include <string>
#include <mpi.h>
#include <logger.hpp>
#include <types.hpp>
#include <graph_types.hpp>
#include <dfs_code.hpp>
#include <graph_output.hpp>
#include <utils.hpp>
#include <deque>
#include <vector>
#include <map>
#include <set>
#include <omp.h>
#include <sys/time.h>

using std::string;

namespace algs {

class global_support_handler {

protected:

  typedef enum { MSG_GLOBAL_SUPPORT = 27, MSG_SUPPORT_DATA = 28} MSG_TAG;

  inline std::string get_msg_tag(MSG_TAG tg) {
    switch(tg) {
    case MSG_GLOBAL_SUPPORT:
      return "MSG_GLOBAL_SUPPORT";
    case MSG_SUPPORT_DATA:
      return "MSG_SUPPORT_DATA";
    default:
      std::stringstream ss;
      ss << "UNKNOWN(" << tg << ")" << endl;
      return ss.str();
    }
  } // get_msg_tag

  MPI_Request request;
  Logger *logger;
  int rank, numtasks, numthreads, minsupp, num_partitions;
  int dfscode_size;
  int bcast_root;
  int rank_finish_count;
  int frequent_patterns;
  bool local_supports_sendrecv_done;
  int patterns_received;
  bool global_support_optimization; //, estimate_with_external_maps;
  std::string vmap_filename;

  std::map<string, std::vector<int> > supports_locally_frequent;
  std::map<string, std::vector<int> > supports_locally_not_frequent;
  std::map<string, std::vector<int> > external_supports, internal_supports;
  std::vector<std::map<string, std::vector<int> > > supports_locally_frequent_threads;
  std::vector<std::map<string, std::vector<int> > > supports_locally_not_frequent_threads;
  std::vector<std::map<string, std::vector<int> > > external_supports_threads, internal_supports_threads;

  std::map<string, std::vector<std::set<int> > > vertex_mappings;  //both locally frequent and not frequent
  std::map<string, std::vector<std::set<int> > > external_mappings;  //dfscode-> vertex_id -> global_vertex_ids
  std::map<string, std::vector<std::map<int, std::set<int> > > > external_mappings_per_part;  //dfscode-> vertex_id -> part_id -> global_vertex_ids
  std::map<string, int> vertex_mappings_file_pos;
  std::map<string, int> vertex_mappings_size_in_file;
  std::map<string, std::vector<int> > global_supports; //dfscode string, support of vertices
  std::map<string, std::vector<int> > global_supports_internal; //dfscode string, support of vertices
  std::map<string, std::vector<int> > global_supports_external; //dfscode string, support of vertices
  std::map<string, std::vector<int> > global_supports_external_max; //dfscode string, max support of vertices
  std::map<int, std::vector<types::DFSCode> > hashed_dfscodes;
  std::vector<types::DFSCode> globally_frequent_dfscodes_hashed;
  std::vector<types::DFSCode> globally_frequent_dfscodes; // includes only locally present ones
  std::vector<types::DFSCode> all_globally_estimated_dfscodes; //includes all, whether locally present or not
  std::map<string, int> marked_actual_globally_frequent_patterns; // if deemed to be actually globally frequent, 1 otherwise 0
  std::map<string, int> final_dfscode_supports; //includes all, whether locally present or not
  std::map<string, int> local_supports_globally_frequent; // local support of globally frequent patterns
  std::map<string, std::vector<std::vector<std::set<int> > > > stored_dfscode_external_vertex_mappings;

  double hash_sendrecv_time;
  double compute_support_time;
  double unfinished_sendrecv_time;
  double global_pattern_sendrecv_time;

  omp_lock_t lock_gs;
  graph_output *out;

public:

  global_support_handler() {
    omp_init_lock(&lock_gs);
  }
  void init_global_support_handler(unsigned int rank, unsigned int numtasks, unsigned int numthreads);

  virtual ~global_support_handler() {
    omp_destroy_lock(&lock_gs);
  }

  int get_globally_frequent();
  virtual void set_min_support(int minsupp);
  void set_num_partitions(int num_part);
  void set_optimization_flags(bool global_support_optimization);
  virtual void set_graph_output(graph_output* gout);
  void set_vmap_filename(std::string filename);
  void insert_local_support_data(types::DFSCode DFS_CODE, std::vector<int> vertex_supp, bool is_locally_frequent);
  void insert_local_support_data_no_lock(types::DFSCode DFS_CODE, std::vector<int> vertex_supp, bool is_locally_frequent);
  void remove_local_support_data(types::DFSCode DFS_CODE, bool is_locally_frequent);
  void insert_local_support_for_globally_frequent(types::DFSCode DFS_CODE, int sup);
  int get_local_support_for_globally_frequent(types::DFSCode DFS_CODE);

  int compute_local_support_and_store_mappings(types::Embeddings &embeddings, types::DFSCode DFS_CODE, bool store_in_file, std::vector<int> &vertex_support);
  int compute_local_support(types::Embeddings &embeddings, types::DFSCode DFS_CODE, std::vector<int> &vertex_support);
  int compute_local_support_and_store_external_mappings(types::Graph &graph, types::Embeddings &embeddings, types::DFSCode DFS_CODE, std::vector<int> &vertex_support);
  int compute_local_support_and_store_int_ext_supports(types::Graph &graph, types::Embeddings &embeddings, types::DFSCode DFS_CODE, std::vector<int> &vertex_support);
  int compute_local_support_and_store_int_ext_supports_no_lock(types::Graph &graph, types::Embeddings &embeddings, types::DFSCode DFS_CODE, std::vector<int> &vertex_support);
  int compute_local_support_and_store_mappings_int_ext_supports(types::Graph &graph, types::Embeddings &embeddings, types::DFSCode DFS_CODE, bool store_in_file, std::vector<int> &vertex_support);

  //void send_buffered_dfscode_supports(std::vector<types::DFSCode> &dfscodes_vec, int dest);
  void gather_and_compute_global_supports(int dfscode_size);
  void gather_estimate_compute_global_supports(int dfscode_size, types::Graph &graph);
  void gather_and_compute_global_supports_with_mappings(int dfscode_size, types::Graph &graph);
  void gather_locally_frequent_patterns_with_no_estimation(int dfscode_size);
  void hash_and_sendrecv_locally_frequent_patterns();
  void hash_and_sendrecv_locally_frequent_patterns_with_int_ext_map_cardinalities();
  void compute_globally_frequent_patterns();
  void compute_globally_frequent_patterns_with_int_ext_map_cardinalities();
  void compute_globally_frequent_patterns_with_mappings(types::Graph &graph);
  int process_non_frequent_local_support_requests(int *buffer, int len, int num_dfscodes, int dest);
  int process_non_frequent_local_support_requests(int *buffer, int len, int num_dfscodes, int dest, int *data, int &datasize);
  int process_non_frequent_local_support_requests_with_int_ext_map_cardinalities(int *buffer, int len, int num_dfscodes, int dest, int *data, int &datasize);
  void process_received_data(int *&recv_buf, int size, int source);
  void process_received_data2(int *recv_buf, int size, int source);
  void process_received_data2_with_int_ext_map_cardinalities(int *recv_buf, int size, int source);
  void sendrecv_unfinished_frequent_patterns();
  void sendrecv_unfinished_frequent_patterns_with_int_ext_map_cardinalities();
  void sendrecv_unfinished_frequent_patterns_with_int_ext_map_cardinalities2();
  void send_globally_frequent_patterns();
  void send_globally_frequent_and_estimated_patterns();
  bool remove_globally_frequent_dfscode(types::DFSCode dfscode);
  void send_globally_frequent_patterns_to_rank_group();
  void process_globally_frequent_patterns(int *buffer, int len, int num_dfscodes, int dest);
  void process_globally_frequent_patterns_with_frequent_flags(int *buffer, int len, int num_dfscodes, int dest);
  void process_globally_frequent_patterns2(int *buffer, int len, int num_dfscodes);
  void get_all_locally_frequent_patterns(std::vector<types::DFSCode> &all_locally_frequent);
  void get_all_locally_not_frequent_patterns(std::vector<types::DFSCode> &all_locally_not_frequent);
  void get_all_locally_frequent_patterns_with_internal_supports(std::vector<types::DFSCode> &all_locally_frequent);
  void get_remaining_frequent_patterns(std::vector<types::DFSCode> &remaining_estimated);
  void get_all_globally_estimated_patterns(std::vector<types::DFSCode> &globally_frequent);
  void get_all_globally_estimated_patterns(std::vector<types::DFSCode> &globally_frequent, int max_support);
  void remove_all_globally_estimated_patterns();
  void set_dfscode_support(types::DFSCode dfscode, int support);
  void populate_globally_frequent_dfscodes(std::vector<types::DFSCode> &dfscodes_to_process);
  int get_total_frequent_dfscodes();
  void combine_threads_map_vectors();
  bool get_stored_vertex_mappings(types::DFSCode, bool, std::vector<std::set<int> >&);
  void get_stored_dfscode_external_vertex_mappings(std::map<string, std::vector<std::vector<std::set<int> > > > &vertex_mappings);
  void clear_vertex_mappings();
  void clear_stored_dfscode_external_vertex_mappings();

  void forward_token();
  //void send_fin_to_all();

  void print_profile_info();
  void print_globally_frequent_dfscodes();
  void print_average_vertex_mappings_info();

  void bcast_msg(int *buffer, int length, int root);
  void all_gather_msg(int *sendbuf, int *recvbuf, int length);
  void all_gatherv_msg(int *sendbuf, int *recvbuf, int sendlen, int *recvlen);
  void all_gather_msg(int *sendbuf, int *recvbuf, int length, MPI_Comm &rank_group);
  void all_gatherv_msg(int *sendbuf, int *recvbuf, int sendlen, int *recvlen, MPI_Comm &rank_group);
  void alltoall_msg(int *sendbuf, int *recvbuf, int length);
  void alltoallv_msg(int *sendbuf, int *recvbuf, int *sendlen, int *recvlen, int *sdispls, int *rdispls);
  void alltoallv_msg(int *sendbuf, int *recvbuf, int *sendlen, int *recvlen);

  void send_data_msg(int *buffer, int length, int dest_proc);
  void recv_data_msg(int *buffer, int length, int src_proc, MPI_Status &stat);

};

}

#endif
