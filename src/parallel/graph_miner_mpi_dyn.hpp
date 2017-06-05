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
#include <dynamic_load_balancing.hpp>
#include <logger.hpp>
#include <deque>

namespace GRAPH_MINER {

class graph_miner_mpi_dyn : public graph_miner, protected algs::dynamic_load_balancing {
protected:

  int numtasks, rank;
  std::vector<std::deque<types::DFS> > dfs_task_queue;       //keep the sibling extensions for each level
  int task_split_level, task_split_threshold;
  bool is_working;
  int embeddings_regeneration_level;
  int current_dfs_level;       //keep track of the level of the candidate tree

  //graph_miner methods
  void project(Projected &, int dfs_level);
  void regenerate_embeddings(Projected &projected, int dfs_level);
  virtual void run_intern();

  //dynamic load balancing methods
  virtual bool working();
  virtual bool can_split_work();
  virtual void split_work(int* &buffer, int &length);
  virtual void process_received_data(int* buffer, int size);

public:
  graph_miner_mpi_dyn();
  void set_task_split_threshold(int thr);
  virtual void set_load_balance_interval(int i);

};

} // namespace GRAPH_MINER
