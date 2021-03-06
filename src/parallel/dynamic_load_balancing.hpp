/*
 * This source is part of the single graph mining algorithm.
 *
 * Copyright 2013 Robert Kessl
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

#ifndef __DYNAMIC_LOAD_BALANCING_HPP__
#define __DYNAMIC_LOAD_BALANCING_HPP__

#include <string>
#include <mpi.h>
#include <logger.hpp>
#include <types.hpp>
#include <graph_types.hpp>
#include <dfs_code.hpp>
#include <deque>

using std::string;

namespace algs {

class dynamic_load_balancing {

protected:

  typedef enum { RT_TOKEN = 0, RT_WORK_REQUEST = 1, RT_WORK_RESPONSE = 2, RT_COMPUTATION_END = 3, RT_DATA = 4 } REQUEST_TYPE;
  typedef enum { BLACK = 5, WHITE = 6 } TOKEN_COLOR;
  typedef enum { MSG_DATA = 10, MSG_REQUEST = 11} MSG_TAG; //, TAG_TOKEN=12
  typedef enum { ARR = 12, RP = 13} DYN_SCHEME;

  inline std::string get_msg_tag(MSG_TAG tg) {
    switch(tg) {
    case MSG_DATA:
      return "MSG_DATA";
    case MSG_REQUEST:
      return "MSG_REQUEST";
    //case TAG_TOKEN:
    //return "TAG_TOKEN";
    default:
      std::stringstream ss;
      ss << "UNKNOWN(" << tg << ")" << endl;
      return ss.str();
    }
  } // get_msg_tag

  inline string get_token_color(TOKEN_COLOR tcol) {
    switch(tcol) {
    case WHITE:
      return "WHITE";
    case BLACK:
      return "BLACK";
    default:
      std::stringstream ss;
      ss << "UNKNOWN(" << tcol << ")" << endl;
      return ss.str();
    }
  } // get_token_color


  inline string get_request_type(REQUEST_TYPE rt) {
    switch(rt) {
    case RT_TOKEN:
      return "RT_TOKEN";
    case RT_WORK_REQUEST:
      return "RT_WORK_REQUEST";
    case RT_WORK_RESPONSE:
      return "RT_WORK_RESPONSE";
    case RT_COMPUTATION_END:
      return "RT_COMPUTATION_END";
    case RT_DATA:
      return "RT_DATA";
    default:
      std::stringstream ss;
      ss << "UNKNOWN(" << rt << ")" << endl;
      return ss.str();
    }
  }

  Logger *logger;
  int rank, numtasks;
  //std::vector<std::deque<types::DFS> >* dfs_task_queue; //task queue reference

  TOKEN_COLOR my_color;
  bool has_token;
  TOKEN_COLOR token_color;
  int next_work_request;
  bool computation_end;
  bool first_round;
  int scheme;

  unsigned long work_counter;
  int load_balance_interval;
  unsigned long no_work_received_counter;


  MPI::Intracomm communication_object;
  std::vector<int> requested_work_from;

public:
  //dynamic_load_balancing(){}
  void init_dynamic_load_balancing(unsigned int rank, unsigned int numtasks, MPI::Intracomm comm_object);
  bool set_load_balancing_scheme(DYN_SCHEME scheme);

  void send_work_request();
  int check_request(bool blocking);
  void process_request(int source);
  void process_work_split_request(int source);

  bool receive_data(int source, int size);
  void do_dynamic_load_balancing(bool blocking);
  void send_other_msg(int *buffer, int length, int dest_proc);
  void send_data_msg(int *buffer, int length, int dest_proc);
  void recv_other_msg(int *buffer, int length, int src_proc, MPI::Status &stat);
  void recv_data_msg(int *buffer, int length, int src_proc, MPI::Status &stat);

  void process_token();
  void forward_token(TOKEN_COLOR tcol);
  bool dijkstra_is_computation_finished();
  void initiate_termination();
  void load_balance();

  virtual void set_load_balance_interval(int i);

  //pure virtual functions, to be implemented by the subclass
  virtual bool can_split_work() = 0;
  virtual void split_work( int* &buffer, int &length) = 0;
  virtual bool working() = 0;
  virtual void process_received_data(int* buffer, int size) = 0;

};

}

#endif
