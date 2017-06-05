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

#include <graph_miner.hpp>
#include <dynamic_hybrid_load_balancing_priv_params.hpp>
#include <dynamic_hybrid_threads_load_balancing_priv_params.hpp>
#include <logger.hpp>
#include <deque>

using algs::Thread_private_data;

namespace GRAPH_MINER {

class graph_miner_mpi_omp_hybrid : public graph_miner, protected algs::dynamic_load_balancing, protected algs::dynamic_threads_load_balancing {

protected:

  int numtasks, rank;
  int threads_per_rank;
  int task_split_threshold;
  std::vector<int> task_split_level;
  bool is_working;
  int donor_thread;
  std::vector<bool> thread_is_working;
  std::vector<int> frequent_patterns_count;

  std::vector<types::DFSCode> DFS_CODE_V;
  std::vector<DFSCode> DFS_CODE_IS_MIN_V;
  std::vector<Graph> GRAPH_IS_MIN_V;

  std::vector<std::vector<std::deque<types::DFS> > > dfs_task_queue;       //keep the sibling extensions for each level and for each thread
  std::vector<std::vector<std::deque<types::DFS> > > dfs_task_queue_shared;       //keep a separate queue for sharing work
  std::vector<std::vector<std::vector<std::deque<types::DFS> > > > global_shared_queue; //shared queue for global work, one for each rank

  std::vector<int> embeddings_regeneration_level;       //for each thread
  std::vector<int> current_dfs_level;                   //keep track of the level of the candidate tree for each thread

  //graph_miner methods
  bool is_min(int thread_id);
  bool project_is_min(int thread_id, types::Projected &);

  bool is_min(Thread_private_data &);
  bool project_is_min(Thread_private_data &, types::Projected &);

  void project(Projected &, int dfs_level, Thread_private_data &gprv);
  void regenerate_embeddings(Projected &projected, int dfs_level, Thread_private_data &gprv);
  void print_threads_status(int thread_id);
  virtual unsigned int support(Projected &projected, Thread_private_data &gprv);
  virtual void run_intern();

  //dynamic load balancing methods
  virtual bool working();
  virtual void initiate_global_split_work(int requester_id, Thread_private_data &gprv);
  virtual void complete_global_split_work(int requester_id, Thread_private_data &gprv); //process load balancing
  virtual void process_received_data(int* buffer, int size, Thread_private_data &gprv);

  virtual bool all_threads_idle();
  virtual bool thread_working(Thread_private_data &gprv);
  virtual bool can_thread_split_work(Thread_private_data &gprv);
  virtual void thread_split_work(int requesting_thread, int &length, Thread_private_data &gprv);
  virtual void thread_process_received_data(Thread_private_data &gprv);
  virtual void thread_start_working();

  virtual void thread_split_global_work(int requester_rank_id, Thread_private_data &gprv);
  virtual void complete_global_work_split_request(int requester_rank_id, Thread_private_data &gprv); //threads load balancing

public:
  graph_miner_mpi_omp_hybrid(int threads_per_rank);
  virtual ~graph_miner_mpi_omp_hybrid();
  void set_task_split_threshold(int thr);
  virtual void set_load_balance_interval(int i);
  void set_num_threads(int t);
  int get_frequent_patterns_count();
  void print_threads_frequent_patterns();
};


} // namespace GRAPH_MINER
