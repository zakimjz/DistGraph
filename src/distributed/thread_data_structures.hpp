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

#ifndef __THREAD_DATA_STRUCTURES_HPP__
#define __THREAD_DATA_STRUCTURES_HPP__

#include <types.hpp>
#include <graph_types.hpp>
#include <dfs_code.hpp>
#include <vector>
#include <deque>

namespace algs {

struct Thread_private_data {

  int thread_id;
  int task_split_level, embeddings_regeneration_level, current_dfs_level;

  types::DFSCode DFS_CODE;
  types::DFSCode DFS_CODE_IS_MIN;
  types::Graph GRAPH_IS_MIN;
  std::vector<std::deque<types::DFS> >  dfs_task_queue;
  std::deque<types::DFSCode> dfscodes_to_process;

};

}

#endif
