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

#include <graph_miner2.hpp>
#include <thread_data_structures.hpp>
#include <external_neighbor_handler.hpp>
#include <global_support_handler.hpp>
#include <logger.hpp>
#include <deque>
#include <fstream>

using algs::Thread_private_data;

namespace GRAPH_MINER {

class graph_miner_mpi_omp_part : public graph_miner, protected algs::external_neighbor_handler {

protected:

  int numtasks, rank;
  int threads_per_rank, max_threads_per_rank;
  int num_partitions;
  std::vector<int> frequent_patterns_count;
  omp_lock_t lock;
  int max_edge_size;
  int dfs_level_cutoff;
  int single_edge_dfscodes;
  std::string vmap_filename_prefix;

  int num_false_positives;
  bool has_external_neighbor;
  bool purge_externals;
  bool store_embeddings, store_in_file; //or in memory
  bool enable_lower_bound;
  bool no_estimation; //raw computation, all locally frequent patterns are passed thu' false pos test
  std::ofstream time_logger;

  Embeddings_map3 first_embeddings;
  std::vector<types::DFSCode> dfscodes_to_process;
  std::vector<std::vector<types::DFSCode> > dfscodes_to_process_threads_shared;

  std::vector<types::DFSCode> DFS_CODE_V;
  std::vector<DFSCode> DFS_CODE_IS_MIN_V;
  std::vector<Graph> GRAPH_IS_MIN_V;

  std::vector<std::vector<std::deque<types::DFS> > > dfs_task_queue;       //keep the sibling extensions for each level and for each thread
  std::vector<std::vector<std::deque<types::DFS> > > dfs_task_queue_shared;       //keep a separate queue for sharing work

  std::vector<int> embeddings_regeneration_level;       //for each thread
  std::vector<int> current_dfs_level;       //keep track of the level of the candidate tree for each thread

  algs::global_support_handler gs_handler;

  //graph_miner methods
  bool is_min(int thread_id);
  bool project_is_min(int thread_id, types::Projected &);

  bool is_min(Thread_private_data &);
  bool project_is_min(Thread_private_data &, types::Projected &);
  void report(types::DFSCode &code, unsigned int sup);

  int create_first_embeddings();
  void populate_first_dfscodes();
  void remove_non_frequent_first_embeddings();
  void populate_globally_frequent_dfscodes();
  bool is_frequent_single_edge(types::Edge e);
  void get_all_embeddings(Edge first_edge, std::vector<types::Embeddings2> &embeddings,  types::DFSCode DFS_CODE);
  void expand(Embeddings &, int dfs_level, Thread_private_data &gprv);
  void regenerate_embeddings(Embeddings &, int dfs_level, Thread_private_data &gprv);
  void get_vertex_maps(types::DFSCode dfscode, std::vector<std::set<int> > &vertex_maps);
  virtual unsigned int max_vertex_support(Embeddings &, Thread_private_data &gprv, std::vector<int> &vertex_support);
  void run_intern();

  //external neighbor handler methods
  virtual int get_num_neighbors(int source, int global_vid);
  virtual int get_num_neighbors(int src, int global_vid, std::set<int> &exclusions);
  virtual int prepare_neighbor_data_to_send(int source, int global_vid, int *buffer, int size);
  virtual void process_received_neighbor_data(int *buffer, int size);
  virtual void process_received_neighbor_data(int *buffer, int size, int num_neighbors);
  virtual void remove_non_used_pseudo_locals();

  //others
  virtual void compute_exact_support();
  void determine_has_external_neighbor();

public:
  graph_miner_mpi_omp_part(int threads_per_rank, int num_part);
  virtual ~graph_miner_mpi_omp_part();
  void set_time_logger(char *filename);
  virtual void set_min_support(int minsup);
  virtual void set_graph_output(graph_output * gout);
  void set_num_threads(int t);
  void set_max_edge_size(int size);
  void set_vmap_filename_prefix(std::string prefix);
  int get_frequent_patterns_count();
  void print_threads_frequent_patterns();

  virtual void run();
};


} // namespace GRAPH_MINER
