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


#ifndef __EXTERNAL_NEIGHBOR_HANDLER_HPP__
#define __EXTERNAL_NEIGHBOR_HANDLER_HPP__

#include <string>
#include <mpi.h>
#include <logger.hpp>
#include <types.hpp>
#include <graph_types.hpp>
#include <dfs_code.hpp>
#include <utils.hpp>
#include <deque>
#include <set>
#include <map>
#include <omp.h>

using std::string;

namespace algs {

class external_neighbor_handler {

protected:

  typedef enum {RT_NEIGH_REQ = 45, RT_NEIGH_DATA = 46, RT_NEIGH_FIN = 47, RT_NEIGH_TOKEN = 48} REQUEST_TYPE;
  typedef enum { MSG_NEIGH_DATA = 20} MSG_TAG;

  inline std::string get_msg_tag(MSG_TAG tg) {
    switch(tg) {
    case MSG_NEIGH_DATA:
      return "MSG_NEIGH_DATA";
    default:
      std::stringstream ss;
      ss << "UNKNOWN(" << tg << ")" << endl;
      return ss.str();
    }
  } // get_msg_tag

  inline string get_request_type(REQUEST_TYPE rt) {
    switch(rt) {
    case RT_NEIGH_REQ:
      return "RT_NEIGH_REQ";
    case RT_NEIGH_DATA:
      return "RT_NEIGH_DATA";
    case RT_NEIGH_TOKEN:
      return "RT_NEIGH_TOKEN";
    case RT_NEIGH_FIN:
      return "RT_NEIGH_FIN";
    default:
      std::stringstream ss;
      ss << "UNKNOWN(" << rt << ")" << endl;
      return ss.str();
    }
  }

  MPI_Request request;

  Logger *logger;
  int rank, numtasks, numthreads, num_partitions;
  int neighbors_fetched;
  std::map<int, int> external_neighbor_requests;
  std::map<std::string, std::set<int> > requests_for_dfscodes;
  std::vector<std::map<std::string, std::set<int> > > requests_for_dfscodes_threads;
  std::map<int, std::set<int> > pseudo_locals_on_extensions_per_level;       //dfs_level-> global ids// Graph:: vertex_part_id != orig_part_id
  std::set<int> pseudo_locals_on_extensions;       // Graph:: vertex_part_id != orig_part_id
  std::map<int, bool> is_external_neighbor_fetched;
  std::map<int, std::set<int> > processed_external_neighbors;
  omp_lock_t lock_nr;
  std::vector<bool> threads_waiting;
  int ranks_finished;
  bool has_token;
  bool computation_finished;

public:

  external_neighbor_handler(){
    omp_init_lock(&lock_nr);
  }

  virtual ~external_neighbor_handler(){
    omp_destroy_lock(&lock_nr);
  }

  void init_external_neighbor_handler(unsigned int rank, unsigned int numtasks, unsigned int numthreads);
  void set_num_partitions(int num_part);

  void add_neighbor_request(int part_id, int global_id);
  void add_neighbor_request(types::DFSCode dfscode, int part_id, int global_id);
  void add_neighbor_request_no_lock(types::DFSCode dfscode, int part_id, int global_id);
  void remove_neighbor_requests_for_non_frequent_dfscodes(std::vector<types::DFSCode> &dfscodes);
  void send_neighbor_requests();
  void send_neighbor_request(int global_vid, int dest);
  void insert_pseudo_local(int local_vid);
  void insert_pseudo_local_per_level(int level, int local_vid);
  void ext_handler_combine_threads_map_vectors();

  void union_neighbor_requests_from_rank_group();
  void union_used_pseudo_locals_in_rank_group(std::set<int>& used_pseudo_locals);

  int check_message(int &size);
  void handle_external_neighbors();
  void handle_external_neighbors(int requester);
  void send_recv_external_neighbors();

  void process_request(int source, int size);
  void process_neighbor_request(int source, int global_vid);
  void process_neighbor_request(int source, int *req, int len);
  bool receive_neighbor_data(int source, int *data, int size);

  void remove_neighbor_request(int global_id);
  bool is_fetched(int global_id);
  bool is_fetched(std::set<int> global_ids);
  bool all_fetched();
  void print_neighbor_request_list();

  void forward_token();
  void send_fin_to_all();

  void send_msg(int *buffer, int length, int dest_proc);
  void recv_msg(int *buffer, int length, int src_proc, MPI_Status &stat);
  void all_gather_msg(int *sendbuf, int *recvbuf, int length, MPI_Comm rank_group);
  void all_gatherv_msg(int *sendbuf, int *recvbuf, int sendlen, int *recvlen, MPI_Comm rank_group);
  void alltoall_msg(int *sendbuf, int *recvbuf, int length, MPI_Comm rank_group);
  void alltoallv_msg(int *sendbuf, int *recvbuf, int *sendlen, int *recvlen, MPI_Comm rank_group);

  //pure virtual functions, to be implemented by the subclass
  virtual int get_num_neighbors(int source, int global_vid) = 0;
  virtual int get_num_neighbors(int source, int global_vid, std::set<int> &exclusions) = 0;
  virtual int prepare_neighbor_data_to_send(int source, int global_vid, int* buffer, int size) = 0;
  virtual void process_received_neighbor_data(int* buffer, int size) = 0;
  virtual void process_received_neighbor_data(int* buffer, int size, int num_neighbors) = 0;
  virtual void remove_non_used_pseudo_locals() = 0;
};

}

#endif
