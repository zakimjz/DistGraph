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


#ifndef __DYNAMIC_THREADS_LOAD_BALANCING_HPP__
#define __DYNAMIC_THREADS_LOAD_BALANCING_HPP__

#include <string>
#include <mpi.h>
#include <logger.hpp>
#include <types.hpp>
#include <graph_types.hpp>
#include <dfs_code.hpp>
#include <deque>
#include <vector>
#include <omp.h>
#include <thread_data_structures.hpp>

using std::string;

namespace algs {

class dynamic_threads_load_balancing {

protected:

  typedef enum { RT_WORK_REQUEST = 0, RT_WORK_RESPONSE = 1, GLOBAL_WORK_REQUEST = 2, GLOBAL_WORK_RESPONSE = 3} REQUEST_TYPE;
  typedef enum { MSG_DATA = 5, MSG_REQUEST = 6} MSG_TAG; //, TAG_TOKEN=12
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


  inline string get_request_type(REQUEST_TYPE rt) {
    switch(rt) {
    case RT_WORK_REQUEST:
      return "RT_WORK_REQUEST";
    case RT_WORK_RESPONSE:
      return "RT_WORK_RESPONSE";
    case GLOBAL_WORK_REQUEST:
      return "GLOBAL_WORK_REQUEST";
    case GLOBAL_WORK_RESPONSE:
      return "GLOBAL_WORK_RESPONSE";
    default:
      std::stringstream ss;
      ss << "UNKNOWN(" << rt << ")" << endl;
      return ss.str();
    }

  }

  Logger *logger;
  int rank, num_threads;
  std::vector<int> next_work_request;
  int scheme;

  std::vector<std::deque<int> > message_queue;
  std::vector<std::vector<int> > requested_work_from;
  std::map<int, int> global_work_response_counter;
  int all_idle_work_request_counter;
  omp_lock_t lock;
  std::vector<omp_lock_t*> qlock;

  void init_dynamic_threads_load_balancing(int rank, int num_threads);
  void reinitialize_dynamic_threads_load_balancing();
  bool set_load_balancing_scheme(DYN_SCHEME scheme);
  void set_num_threads(int num_threads);
  bool all_idle_work_requests_collected();
  void reset_idle_work_requests_counter();

  void send_work_request(Thread_private_data &gprv);
  void send_work_request(Thread_private_data &gprv, int dest);
  int check_request(Thread_private_data &gprv);
  void process_request(int source, Thread_private_data &gprv);
  void process_work_split_request(int source, Thread_private_data &gprv);
  bool receive_data(int source, int size, Thread_private_data &gprv);
  void send_msg(int *buffer, int length, int src_thr, int dest_thr);
  void recv_msg(int *buffer, int length, int thr, int originating_thr);

  //hybrid load balancing functions
  void send_global_split_requests(int requester_rank_id, Thread_private_data &gprv);
  void process_global_work_split_request(int requester_rank_id, Thread_private_data &gprv);
  void process_global_work_split_response(int requester_rank_id, Thread_private_data &gprv);
  void send_work_response_to_all_threads(Thread_private_data &gprv, std::vector<int> thread_has_work);

public:

  dynamic_threads_load_balancing(){
    omp_init_lock(&lock);
  };

  virtual ~dynamic_threads_load_balancing(){
    omp_destroy_lock(&lock);
  }

  void threads_load_balance(Thread_private_data &gprv);

  //pure virtual functions, to be implemented by the subclass
  virtual bool can_thread_split_work(Thread_private_data &gprv) = 0;
  virtual void thread_split_work(int requesting_thread, int &length, Thread_private_data &gprv) = 0;
  virtual void thread_process_received_data(Thread_private_data &gprv) = 0;
  virtual bool thread_working(Thread_private_data &gprv) = 0;
  virtual bool all_threads_idle() = 0;
  virtual void thread_start_working() = 0;

  //hybrid load balancing functions
  virtual void thread_split_global_work(int requester_rank_id, Thread_private_data &gprv) = 0;
  virtual void complete_global_work_split_request(int requester_rank_id, Thread_private_data &gprv) = 0;
};

}

#endif
