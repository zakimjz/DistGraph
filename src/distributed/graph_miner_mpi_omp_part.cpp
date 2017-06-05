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

#include <graph_miner_mpi_omp_part.hpp>
#include <dbio.hpp>
#include <iterator>
#include <stdlib.h>
#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <math.h>
#include <omp.h>
#include <cassert>

using namespace std;

namespace GRAPH_MINER {

graph_miner_mpi_omp_part::graph_miner_mpi_omp_part(int th_per_rank, int num_part) : graph_miner()
{

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  max_edge_size = 100000;   //set a large number for max pattern size, such that we don't expect to find patterns of that size
  dfs_level_cutoff = 100000; //dfs level in the pattern lattice corresponds to the pattern size, should also be set to a large number
  single_edge_dfscodes = 0;
  num_partitions = num_part;
  DEBUG(*(graph_miner::logger), " omp_max_threads = " << omp_get_max_threads());
  if(th_per_rank < omp_get_max_threads()) {
    max_threads_per_rank = omp_get_max_threads();
    threads_per_rank = max_threads_per_rank;
  }else{
    max_threads_per_rank = th_per_rank;
    threads_per_rank = th_per_rank;
  }
  num_false_positives = 0;

  store_embeddings = false; //set to true for both memory and file based storage
  store_in_file = false;   // true = file, false = memory --- based storage
  enable_lower_bound = true; //enable lower bound estimation of support, upper bound is already enabled by default

  no_estimation = false;  //disable both lower and upper bound support estimation

  has_external_neighbor = false; //if the partitions have external edges (recomputed later from the input), so setting it to true won't have any affect
  purge_externals = false;  // enable external neighbor purging/shrinking of local partition after each iteration

  for(int i = 0; i<threads_per_rank; i++) {

    frequent_patterns_count.push_back(0);
    embeddings_regeneration_level.push_back(0);
    current_dfs_level.push_back(0);

    types::DFSCode dfscode;
    DFS_CODE_V.push_back(dfscode);
    DFS_CODE_IS_MIN_V.push_back(dfscode);

    types::Graph gr;
    GRAPH_IS_MIN_V.push_back(gr);

    std::vector<std::deque<types::DFS> > tmp;
    dfs_task_queue.push_back(tmp);
    dfs_task_queue_shared.push_back(tmp);

    std::vector<types::DFSCode> tmp2;
    dfscodes_to_process_threads_shared.push_back(tmp2);

  }

  omp_init_lock(&lock);

  //char time_logfile[100];
  //sprintf(time_logfile, "time-log-p%d", rank);
  //time_logger.open(time_logfile);

  init_external_neighbor_handler(rank, numtasks, threads_per_rank);
  external_neighbor_handler::set_num_partitions(num_part);
  gs_handler.init_global_support_handler(rank, numtasks, threads_per_rank);
  gs_handler.set_num_partitions(num_part);
  gs_handler.set_optimization_flags(enable_lower_bound);
}

graph_miner_mpi_omp_part::~graph_miner_mpi_omp_part(){

  omp_destroy_lock(&lock);
  if(time_logger.is_open())
    time_logger.close();

}

void graph_miner_mpi_omp_part::set_time_logger(char *filename){
  if(time_logger.is_open())
    time_logger.close();
  time_logger.open(filename);
}

void graph_miner_mpi_omp_part::set_num_threads(int t){
  omp_set_num_threads(t);
}

void graph_miner_mpi_omp_part::set_max_edge_size(int size){
  max_edge_size = size;
}

void graph_miner_mpi_omp_part::set_min_support(int minsup){
  minimal_support = minsup;
  gs_handler.set_min_support(minsup);
}

void graph_miner_mpi_omp_part::set_graph_output(graph_output * gout)
{
  output = gout;
  if(gout)
    gs_handler.set_graph_output(gout);
}

int graph_miner_mpi_omp_part::get_frequent_patterns_count(){
  return gs_handler.get_globally_frequent();
}

void graph_miner_mpi_omp_part::set_vmap_filename_prefix(std::string prefix){
  vmap_filename_prefix = prefix;
  gs_handler.set_vmap_filename(prefix);
}

void graph_miner_mpi_omp_part::print_threads_frequent_patterns(){
  stringstream ss;
  for(int i = 0; i < threads_per_rank; i++)
    ss << " tid = " << i << " patterns = " << frequent_patterns_count[i];
  INFO(*(graph_miner::logger), ss.str());
}

void graph_miner_mpi_omp_part::report(types::DFSCode &code, unsigned int sup)
{
  //omp_set_lock(&lock);
  output->output_graph(code, sup);
  //omp_unset_lock(&lock);
}

//support function for a single large graph, computes the maximum count of a node in the embeddings
unsigned int
graph_miner_mpi_omp_part::max_vertex_support(Embeddings &embeddings, Thread_private_data &gprv, std::vector<int> &vertex_support)
{
  std::map<unsigned int, map<unsigned int, unsigned int> > node_id_counts;

  int maxtoc = 0;
  int thread_id = gprv.thread_id; //omp_get_thread_num();
  //iterated through the all the embeddings
  for(Embeddings::iterator cur = embeddings.begin(); cur != embeddings.end(); ++cur) {

    Emb *em = &(*cur);
    int dfsindex = gprv.DFS_CODE.size() - 1;
    while(em != NULL) {
      if(gprv.DFS_CODE[dfsindex].to > gprv.DFS_CODE[dfsindex].from) {    //forward edge
        node_id_counts[gprv.DFS_CODE[dfsindex].to][(*em).edge.to]++;
      }
      if(!em->prev) {
        node_id_counts[gprv.DFS_CODE[dfsindex].from][(*em).edge.from]++;
      }

      if(maxtoc < gprv.DFS_CODE[dfsindex].to)
        maxtoc = gprv.DFS_CODE[dfsindex].to;

      em = em->prev;
      dfsindex--;
    }
  }

  unsigned int max = 0;
  for(std::map<unsigned int, map<unsigned int, unsigned int> >::iterator it = node_id_counts.begin(); it != node_id_counts.end(); it++) {
    if((it->second).size() > max)
      max = (it->second).size();
  }

  for(int i = 0; i<=maxtoc; i++)
    vertex_support.push_back(node_id_counts[i].size());

  return max;
}


void graph_miner_mpi_omp_part::expand(Embeddings &embeddings, int dfs_level, Thread_private_data &gprv)
{
  // Check if the pattern is frequent enough.
  std::vector<int> vertex_support;
  int thread_id = gprv.thread_id; //omp_get_thread_num();

  //omp_set_lock(&lock);
  //DEBUG(*(graph_miner::logger), "prior to support compute: thread " << thread_id <<" executing expand for code: " << gprv.DFS_CODE.to_string() << "; embeddings size: " << embeddings.size() );
  //omp_unset_lock(&lock);

  //if(gprv.DFS_CODE.to_string() == "F(0 1 1 0 1)")
  //embeddings.print_global_vid(graph);

  //allow extensions upto dfs_level_cutoff without support computation and mindfs checking
  unsigned int sup;
  if(dfs_level >= dfs_level_cutoff) {
    if(is_min(gprv) == false)
      return;

    if(store_embeddings) {
      if(enable_lower_bound)
        sup = gs_handler.compute_local_support_and_store_mappings_int_ext_supports(graph, embeddings, gprv.DFS_CODE, store_in_file, vertex_support);
      else
        sup = gs_handler.compute_local_support_and_store_mappings(embeddings, gprv.DFS_CODE, store_in_file, vertex_support);
    } else{
      if(enable_lower_bound)
        //sup = gs_handler.compute_local_support_and_store_int_ext_supports(graph, embeddings, gprv.DFS_CODE, vertex_support);
        sup = gs_handler.compute_local_support_and_store_int_ext_supports_no_lock(graph, embeddings, gprv.DFS_CODE, vertex_support);
      else
        sup = max_vertex_support(embeddings, gprv, vertex_support);
    }

    DEBUG(*(graph_miner::logger), "DFS Level " << dfs_level << " thread " << thread_id << " executing expand for code: " << gprv.DFS_CODE.to_string() << "; support: " << sup << "; embedddings size = " << embeddings.size() );

    bool is_locally_frequent;
    is_locally_frequent = (sup >= (minimal_support / num_partitions) );
    //gs_handler.insert_local_support_data(gprv.DFS_CODE, vertex_support, is_locally_frequent);
    gs_handler.insert_local_support_data_no_lock(gprv.DFS_CODE, vertex_support, is_locally_frequent);


    DEBUG(*logger, "DFS level = " << dfs_level);
    //omp_set_lock(&lock);
    DEBUG(*(graph_miner::logger), "DFS Level " << dfs_level << " thread " << thread_id << " executing expand for code: " << gprv.DFS_CODE.to_string() << "; support: " << sup << "; embedddings size = " << embeddings.size() );
    //omp_unset_lock(&lock);
  }

  if(gprv.DFS_CODE[dfs_level - 1].is_forward()) {
    for(unsigned int n = 0; n < embeddings.size(); ++n) {
      Edge e = embeddings[n].edge;

      if (graph[e.to].vertex_part_id != rank % num_partitions) {
        //add_neighbor_request(gprv.DFS_CODE, graph[e.to].vertex_part_id, graph[e.to].global_vid);
        add_neighbor_request_no_lock(gprv.DFS_CODE, graph[e.to].vertex_part_id, graph[e.to].global_vid);
      }
      if (graph[e.from].vertex_part_id != rank % num_partitions) {
        //add_neighbor_request(gprv.DFS_CODE, graph[e.from].vertex_part_id, graph[e.from].global_vid);
        add_neighbor_request_no_lock(gprv.DFS_CODE, graph[e.from].vertex_part_id, graph[e.from].global_vid);
      }

      if(purge_externals) {
        if(dfs_level == 1 && graph.is_pseudo_local(e.from)) {
          insert_pseudo_local_per_level(dfs_level, graph[e.from].global_vid);
        }
        if(gprv.DFS_CODE[dfs_level - 1].is_forward() && graph.is_pseudo_local(e.to)) {
          insert_pseudo_local_per_level(dfs_level, graph[e.to].global_vid);
        }
      }

    }
  }

  if (dfs_level == dfs_level_cutoff || dfs_level == max_edge_size)
    return;
  //gprv.frequent_patterns_count++;
  //std::cout << "Frequent DFS Code: "<<gprv.DFS_CODE.to_string() << " sup = "<< sup <<std::endl;

  // In case we have a valid upper bound and our graph already exceeds it,
  // return.  Note: we do not check for equality as the DFS exploration may
  // still add edges within an existing subgraph, without increasing the
  // number of nodes.
  //
  //if(maxpat_max > maxpat_min && DFS_CODE.nodeCount() > maxpat_max) return;

  // We just outputted a frequent subgraph.  As it is frequent enough, so
  // might be its (n+1)-extension-graphs, hence we enumerate them all.

  const RMPath &rmpath = gprv.DFS_CODE.buildRMPath();
  int minlabel = gprv.DFS_CODE[0].fromlabel;
  int maxtoc = gprv.DFS_CODE[rmpath[0]].to; //right-most vertex

  Embeddings_map3 new_fwd_root;
  Embeddings_map2 new_bck_root;
  std::vector<Edge> edges;

  gprv.current_dfs_level = dfs_level;

  for(unsigned int n = 0; n < embeddings.size(); ++n) {

    Emb *cur = &embeddings[n];
    EmbVector history(graph, cur);

    // XXX: do we have to change something here for directed edges?

    // backward
    for(int i = (int)rmpath.size() - 1; i >= 1; --i) {
      Edge e;
      if(get_backward(graph, history[rmpath[i]], history[rmpath[0]], history, e)) {
        if(is_frequent_single_edge(e))
          new_bck_root[gprv.DFS_CODE[rmpath[i]].from][e.elabel].push(e, cur);
      }
    }

    // pure forward
    // FIXME: here we pass a too large e->to (== history[rmpath[0]]->to
    // into get_forward_pure, such that the assertion fails.
    //
    // The problem is:
    // history[rmpath[0]]->to > graph_database[id].size()

    if(get_forward_pure(graph, history[rmpath[0]], minlabel, history, edges)) {
      for(std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); ++it) {
        //Check if the extension is a frequent single dfs code
        if(is_frequent_single_edge(*it))
          new_fwd_root[maxtoc][it->elabel][graph[it->to].label].push(*it, cur);

      }
    }

    // backtracked forward
    for(int i = 0; i < (int)rmpath.size(); ++i) {
      if(get_forward_rmpath(graph, history[rmpath[i]], minlabel, history, edges)) {
        for(std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); ++it) {
          //Check if the extension is a frequent single dfs code
          if(is_frequent_single_edge(*it))
            new_fwd_root[gprv.DFS_CODE[rmpath[i]].from][it->elabel][graph[it->to].label].push(*it, cur);

        }         // for it
      } // if
    } // for i

  } // for n

  std::deque<types::DFS> tmp;

  if(gprv.dfs_task_queue.size() <= dfs_level) {
    gprv.dfs_task_queue.push_back(tmp);
  }

  // Test all extended substructures.
  // backward
  for(Embeddings_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
    for(Embeddings_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {

      DFS dfs(maxtoc, to->first, -1, elabel->first, -1);
      //DFS dfs(maxtoc, to->first, graph[maxtoc].label, elabel->first, graph[to->first].label);
      gprv.dfs_task_queue[dfs_level].push_back(dfs);
    }
  }

  // forward
  for(Embeddings_riterator3 from = new_fwd_root.rbegin();
      from != new_fwd_root.rend(); ++from) {
    for(Embeddings_iterator2 elabel = from->second.begin();
        elabel != from->second.end(); ++elabel) {
      for(Embeddings_iterator1 tolabel = elabel->second.begin();
          tolabel != elabel->second.end(); ++tolabel) {

        DFS dfs(from->first, maxtoc + 1, -1, elabel->first, tolabel->first);
        //DFS dfs(from->first, maxtoc + 1, graph[from->first].label, elabel->first, tolabel->first);
        gprv.dfs_task_queue[dfs_level].push_back(dfs);
      }
    }
  }

  while(gprv.dfs_task_queue[dfs_level].size() > 0) {

    DFS dfs = gprv.dfs_task_queue[dfs_level].front();
    gprv.dfs_task_queue[dfs_level].pop_front();

    DEBUG(*(graph_miner::logger), "popped dfs = " << dfs.to_string() );

    gprv.current_dfs_level = dfs_level;

    gprv.DFS_CODE.push(dfs.from, dfs.to, dfs.fromlabel, dfs.elabel, dfs.tolabel);

    if(dfs.is_backward())
      expand(new_bck_root[dfs.to][dfs.elabel], dfs_level + 1, gprv);      //Embeddings (Emb vector): each entry contains graph id 0, edge pointer, null Emb
    else
      expand(new_fwd_root[dfs.from][dfs.elabel][dfs.tolabel], dfs_level + 1, gprv);      //Embeddings (Emb vector): each entry contains graph id 0, edge pointer, null Emb

    gprv.DFS_CODE.pop();

  }

  return;
}

/*
 * Regnerate the embeddings for the dfs codes obtained from other processors
 */
void graph_miner_mpi_omp_part::regenerate_embeddings(Embeddings &embeddings, int dfs_level, Thread_private_data &gprv)
{
  // We don't need to check if the pattern is frequent or minimal

  //omp_set_lock(&lock);
  DEBUG(*(graph_miner::logger), "DFS level inside regenerate embeddings = " << dfs_level << " DFS Code " << gprv.DFS_CODE.to_string() << " queue size = " << gprv.dfs_task_queue[dfs_level].size() << " embeddings size = " << embeddings.size());
  //omp_unset_lock(&lock);

  int thread_id = gprv.thread_id; //omp_get_thread_num();

  //iterate for all in the task_queue
  for(int i = 0; gprv.dfs_task_queue[dfs_level].size() > 0; i++) {

    gprv.current_dfs_level = dfs_level;

    types::DFS dfs = gprv.dfs_task_queue[dfs_level].front();
    gprv.dfs_task_queue[dfs_level].pop_front();

    gprv.DFS_CODE.push(dfs.from, dfs.to, dfs.fromlabel, dfs.elabel, dfs.tolabel);

    Embeddings new_root;

    for(unsigned int n = 0; n < embeddings.size(); ++n) {

      Emb *cur = &embeddings[n];
      EmbVector history(graph, cur);

      if(dfs.is_backward() ) {
        Edge e;
        if(get_backward(graph, gprv.DFS_CODE, history, e))
          new_root.push(e, cur);
      }else{
        std::vector<Edge> edges;
        if(get_forward(graph, gprv.DFS_CODE, history, edges)) {
          for(std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); ++it) {
            new_root.push(*it, cur);

            if(purge_externals) {
              if(dfs_level == 1 && graph.is_pseudo_local(it->from)) {
                insert_pseudo_local_per_level(dfs_level, graph[it->from].global_vid);
              }
              if(gprv.DFS_CODE[dfs_level - 1].is_forward() && graph.is_pseudo_local(it->to)) {
                insert_pseudo_local_per_level(dfs_level, graph[it->to].global_vid);
              }
            }

          }
        }
      }
    }

    if( gprv.embeddings_regeneration_level > dfs_level ) {
      regenerate_embeddings(new_root, dfs_level + 1, gprv);
    }else{
      //regeneration of embeddings ended
      //now perform regular extensions with expand function
      expand(new_root, dfs_level + 1, gprv);
    }
    gprv.DFS_CODE.pop();
  }

  return;
}

void graph_miner_mpi_omp_part::get_all_embeddings(Edge first_edge, std::vector<types::Embeddings2> &embeddings,  types::DFSCode DFS_CODE) {

  //embeddings will contain all isomorphisms starting from the first edge
  //one to one mapping from DFS_CODE to embeddings
  types::Embeddings2 t;
  embeddings.push_back(t);
  embeddings[0].push(first_edge, -1);

  types::DFSCode DFS_PREFIX;
  DFS_PREFIX.push_back(DFS_CODE[0]);

  //skip the first edge
  for(int i = 1; i < DFS_CODE.size(); ++i) {

    types::DFS dfs = DFS_CODE[i];
    DFS_PREFIX.push_back(dfs);

    embeddings.push_back(t);
    //get all one edge extensions
    for(int n = 0; n < embeddings[i - 1].size(); n++) {
      EmbVector history(graph, embeddings, i - 1, n);
      if(dfs.is_backward()) {
        Edge e;
        if(get_backward(graph, DFS_PREFIX, history, e))
          embeddings[i].push(e, n);
      }else{                    //forward
        std::vector<Edge> edges;
        if(get_forward(graph, DFS_PREFIX, history, edges)) {
          for(std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); ++it) {
            //prune it->to if not in vertex_mappings[dfs.to]
            embeddings[i].push(*it, n);
          }
        }
      }
    }
  }
}


void graph_miner_mpi_omp_part::get_vertex_maps(types::DFSCode DFS_CODE_full, std::vector<std::set<int> > &vertex_maps){

  types::DFS dfs = DFS_CODE_full[0];
  Embeddings& first_embs = first_embeddings[dfs.fromlabel][dfs.elabel][dfs.tolabel];

  const RMPath &rmpath = DFS_CODE_full.buildRMPath();
  int minlabel = DFS_CODE_full[0].fromlabel;
  int maxtoc = DFS_CODE_full[rmpath[0]].to;

  //std::vector<int> vmap(vertex_sets_size);

  DFSCode DFS_CODE2 = DFS_CODE_full;

  //Expand the embeddings and get all candidate extensions
  for(int c = 0; c < first_embs.size(); c++) {
    std::vector<types::Embeddings2> embs2;
    get_all_embeddings(first_embs[c].edge, embs2, DFS_CODE_full);

    for(int n = 0; n < embs2[DFS_CODE_full.size() - 1].size(); n++) {
      EmbVector history(graph, embs2, DFS_CODE_full.size() - 1, n);
      //get the unique vertex mappings
      for(int k = 0; k < history.size(); k++) {
        Edge e = history[k];
        if(DFS_CODE_full[k].is_forward()) {                                   //forward edge
          vertex_maps[DFS_CODE_full[k].to].insert(e.to);
          //vmap[DFS_CODE[k].to] = e.to;
        }
        if(k == 0) {
          vertex_maps[0].insert(e.from);
          //vmap[0] = e.from;
        }
      }
    }
  }
}

bool graph_miner_mpi_omp_part::is_frequent_single_edge(types::Edge e){
  //Check if the extension is a frequent single dfs code
  return ( (first_embeddings.count(graph[e.from].label) != 0 && first_embeddings[graph[e.from].label].count(e.elabel) != 0
            && first_embeddings[graph[e.from].label][e.elabel].count(graph[e.to].label) != 0)
           || ( first_embeddings.count(graph[e.to].label) != 0 && first_embeddings[graph[e.to].label].count(e.elabel) != 0
                && first_embeddings[graph[e.to].label][e.elabel].count(graph[e.from].label) != 0 ) );
}

int graph_miner_mpi_omp_part::create_first_embeddings(){

  DEBUG(*(graph_miner::logger)," Rank " << rank << " graph size " << graph.size() );

  std::vector<Edge> edges;
  //Embeddings_map3 root;
  //int single_edge_dfscodes = 0; // now a property

  for(unsigned int from = 0; from < graph.size(); ++from) {
    if(get_forward_root(graph, graph[from], edges)) {           // get the edge list of the node graph[from] in graph g
      for(std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); ++it) {

        //embeddings with a single edge
        if(graph[from].label <= graph[it->to].label) {                //check DFS code minimality for first edge
          if(first_embeddings.count(graph[from].label) == 0 || first_embeddings[graph[from].label].count(it->elabel) == 0 || first_embeddings[graph[from].label][it->elabel].count(graph[it->to].label) == 0) {
            single_edge_dfscodes++;
            DEBUG(*(graph_miner::logger), "single edge DFS code : (0,1," << graph[from].label << "," << it->elabel << "," << graph[it->to].label << ")" );
          }
          first_embeddings[graph[from].label][it->elabel][graph[it->to].label].push(*it, 0);                                //embeddings (Emb vector) entry: graph id (always 0 for single graph), edge pointer and null Emb
        }

      }          //for
    }           // if
  }           // for from
  //} // for id

  return single_edge_dfscodes;
}


void graph_miner_mpi_omp_part::populate_first_dfscodes(void){

  for(Embeddings_iterator3 fromlabel = first_embeddings.begin(); fromlabel != first_embeddings.end(); ++fromlabel) {
    for(Embeddings_iterator2 elabel = fromlabel->second.begin(); elabel != fromlabel->second.end(); ++elabel) {
      for(Embeddings_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
        DFS dfs(0, 1, fromlabel->first, elabel->first, tolabel->first);
        types::DFSCode dfscode;
        dfscode.push_back(dfs);
        dfscodes_to_process.push_back(dfscode);
      }     // for tolabel
    }     // for elabel
  }

}

void graph_miner_mpi_omp_part::populate_globally_frequent_dfscodes(void){
  dfscodes_to_process.clear();
  gs_handler.populate_globally_frequent_dfscodes(dfscodes_to_process);
}

void graph_miner_mpi_omp_part::remove_non_frequent_first_embeddings(){
  types::Embeddings_map3 temp_embeddings;
  for(int i = 0; i < dfscodes_to_process.size(); i++) {
    DFS dfs = dfscodes_to_process[i][0];             //single edge code
    temp_embeddings[dfs.fromlabel][dfs.elabel][dfs.tolabel] = first_embeddings[dfs.fromlabel][dfs.elabel][dfs.tolabel];
  }

  first_embeddings.clear();
  first_embeddings = temp_embeddings;
}


void graph_miner_mpi_omp_part::run_intern(void)
{
  int num_dfscodes = dfscodes_to_process.size();

  Thread_private_data gprv;

  INFO(*(graph_miner::logger), "dfs level = " << dfs_level_cutoff << " total dfs codes to process  = " << num_dfscodes);

  int num_th = threads_per_rank;
  if (threads_per_rank > num_dfscodes)
    num_th = num_dfscodes;

   #pragma omp parallel num_threads(num_th) private(gprv)
  {

    timeval start_time, stop_time;
    gettimeofday(&start_time, 0);

    int thread_id = omp_get_thread_num();
    gprv.thread_id = thread_id;
    gprv.current_dfs_level = 0;
    gprv.task_split_level = 0;
    gprv.embeddings_regeneration_level = 0;

    Embeddings emb_prv;

    int start_index = 0, end_index = 0;

    int ranks_per_partition = numtasks / num_partitions;
    int dfscodes_per_rank = (int) ceil(num_dfscodes * 1.0 / ranks_per_partition);

    int start_rank_index = (rank / num_partitions) * dfscodes_per_rank;
    int end_rank_index = start_rank_index + dfscodes_per_rank - 1;

    if (end_rank_index > num_dfscodes - 1)
      end_rank_index = num_dfscodes - 1;

    if(start_rank_index <= end_rank_index) {
      int dfscodes_per_thread = (int) ceil((end_rank_index - start_rank_index + 1) * 1.0 / threads_per_rank);
      start_index = start_rank_index + thread_id * dfscodes_per_thread;
      end_index = start_index + dfscodes_per_thread - 1;
      if (end_index > end_rank_index)
        end_index = end_rank_index;
    }else{
      start_index = 1;
      end_index = 0;
    }

    if(start_index <= end_index) { //has some work
      //omp_set_lock(&lock);
      //DEBUG(*(graph_miner::logger), "rank " << rank << " thread id " << thread_id <<  " dfs level = " << dfs_level_cutoff << "start index = " << start_index << " , end index = " << end_index << " total dfs codes to process  = " << num_dfscodes << endl);
      //omp_unset_lock(&lock);

      int initial_max_dfs_level = 0;
      for(std::vector<types::DFSCode>::iterator it = dfscodes_to_process.begin(); it != dfscodes_to_process.end(); ++it) {
        if(it->size() > initial_max_dfs_level)
          initial_max_dfs_level = it->size();
      }

      gprv.dfs_task_queue.clear();
      for(int i = 0; i<initial_max_dfs_level; i++) {
        std::deque<types::DFS> tmp;
        gprv.dfs_task_queue.push_back(tmp);
      }

      int index = 0;
      gprv.dfscodes_to_process.clear();
      for(std::vector<types::DFSCode>::iterator it = dfscodes_to_process.begin(); it != dfscodes_to_process.end(); ++it) {
        if( index >= start_index && index <= end_index ) {
          gprv.dfscodes_to_process.push_back(*it);
        }
        index++;
      }
      DEBUG(*(graph_miner::logger), " rank " << rank << " thread " << thread_id << " total dfscodes " << gprv.dfscodes_to_process.size());

      while( gprv.dfscodes_to_process.size() > 0) {

        types::DFSCode dfscode = gprv.dfscodes_to_process.front();
        gprv.dfscodes_to_process.pop_front();
        DEBUG(*(graph_miner::logger), " gprv dfscodes to process size = " << gprv.dfscodes_to_process.size());
        for(int i = 0; i< dfscode.size(); i++)
          gprv.dfs_task_queue[i].push_back(dfscode[i]);
        gprv.embeddings_regeneration_level = dfs_level_cutoff - 2;

        DFS dfs = gprv.dfs_task_queue[0].front();
        gprv.dfs_task_queue[0].pop_front();

        DEBUG(*(graph_miner::logger), "rank " << rank << " thread " << thread_id << " popped dfs = " << dfs.to_string() << " regen level = " << gprv.embeddings_regeneration_level );
        gprv.DFS_CODE.push(0, 1, dfs.fromlabel, dfs.elabel, dfs.tolabel);

        gprv.current_dfs_level = 1;
        if(gprv.embeddings_regeneration_level < 1)
          expand(first_embeddings[dfs.fromlabel][dfs.elabel][dfs.tolabel], 1, gprv);                                  //Embeddings (Emb vector): each entry contains graph id 0, edge pointer, null Emb
        else
          regenerate_embeddings(first_embeddings[dfs.fromlabel][dfs.elabel][dfs.tolabel], 1, gprv);

        gprv.current_dfs_level = 0;
        gprv.DFS_CODE.pop();

        if(gprv.dfscodes_to_process.size() == 0) {
          gprv.embeddings_regeneration_level = 0;
          embeddings_regeneration_level[thread_id] = 0;
        }
      }             //while

    }

    DEBUG(*(graph_miner::logger), "Done... rank = " << rank << " thread " << thread_id << " computation ended ");

    gettimeofday(&stop_time, 0);

    if(time_logger.is_open()) {
      //omp_set_lock(&lock);
      DEBUG(*(graph_miner::logger), "total_time rank" << rank << " thread " << thread_id << " : " << utils::get_time_diff(start_time, stop_time));
      time_logger << "dfs level " << dfs_level_cutoff << " rank " << rank << " thread " << thread_id << " total_time " << utils::get_time_diff(start_time, stop_time) << std::endl;
      //omp_unset_lock(&lock);
    }

  }  //pragma omp

  MPI_Barrier(MPI_COMM_WORLD);

} // void graph_miner_mpi_omp_part::run_intern_new(void)


void graph_miner_mpi_omp_part::determine_has_external_neighbor(){
  int ext = 0, sum;
  if(graph.has_ext_neighbor)
    ext = 1;

  MPI_Allreduce(&ext, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(sum > 0)
    has_external_neighbor = true;
  else
    has_external_neighbor = false;
}

void graph_miner_mpi_omp_part::run(){

  timeval start_time, stop_time;
  map<string, double> profile_timings;
  string profiles[] = {"mining", "global support computation", "external neighbor handler", "remove false positives", "others"};

  for(int i = 0; i < 5; i++)
    profile_timings[profiles[i]] = 0.0;

  determine_has_external_neighbor();
  if(has_external_neighbor == false) {
    purge_externals = false;
    enable_lower_bound = false;
  }
  gettimeofday(&start_time, 0);

  create_first_embeddings();
  populate_first_dfscodes();

  gettimeofday(&stop_time, 0);
  profile_timings[profiles[0]] += utils::get_time_diff(start_time, stop_time);

  dfs_level_cutoff = 1;
  //max_edge_size = 3;
  while(dfs_level_cutoff <= max_edge_size ) {

    gettimeofday(&start_time, 0);

    DEBUG(*(graph_miner::logger), "rank = " << rank << " dfs level = " << dfs_level_cutoff << " graph size = " << graph.size());

    DEBUG(*(graph_miner::logger), " dfs level = " << dfs_level_cutoff);
    if(purge_externals) {
      std::set<int> t;
      pseudo_locals_on_extensions_per_level[dfs_level_cutoff] = t;
    }

    if(store_embeddings && store_in_file)
      dbio::remove_vmap_file(vmap_filename_prefix.c_str());
    //dbio::remove_vmap_file(rank, 0);

    run_intern();
    gs_handler.combine_threads_map_vectors();
    ext_handler_combine_threads_map_vectors();

    gettimeofday(&stop_time, 0);
    profile_timings[profiles[0]] += utils::get_time_diff(start_time, stop_time);

    gettimeofday(&start_time, 0);

    if(!has_external_neighbor) {
      gs_handler.gather_and_compute_global_supports(dfs_level_cutoff);
    }else{
      if(no_estimation == false)
        gs_handler.gather_estimate_compute_global_supports(dfs_level_cutoff, graph);
      else
        gs_handler.gather_locally_frequent_patterns_with_no_estimation(dfs_level_cutoff);
    }

    gettimeofday(&stop_time, 0);
    profile_timings[profiles[1]] += utils::get_time_diff(start_time, stop_time);

    //DEBUG(*(graph_miner::logger), " dfs level = "<< dfs_level_cutoff << "total frequent codes = "<<gs_handler.get_total_frequent_dfscodes());

    gettimeofday(&start_time, 0);

    if(has_external_neighbor)             // && !estimate_with_external_maps)
      compute_exact_support();
    if(store_embeddings)
      gs_handler.clear_vertex_mappings();
    gettimeofday(&stop_time, 0);
    profile_timings[profiles[3]] += utils::get_time_diff(start_time, stop_time);

    if(gs_handler.get_total_frequent_dfscodes() == 0) break;

    gettimeofday(&start_time, 0);

    gs_handler.print_globally_frequent_dfscodes();
    populate_globally_frequent_dfscodes();
    if(dfs_level_cutoff == 1)
      remove_non_frequent_first_embeddings();

    gettimeofday(&stop_time, 0);
    profile_timings[profiles[4]] += utils::get_time_diff(start_time, stop_time);

    gettimeofday(&start_time, 0);

    if(has_external_neighbor)
      remove_neighbor_requests_for_non_frequent_dfscodes(dfscodes_to_process);

    if(purge_externals) {
      remove_non_used_pseudo_locals();
    }

    if(has_external_neighbor)
      send_recv_external_neighbors();

    gettimeofday(&stop_time, 0);
    profile_timings[profiles[2]] += utils::get_time_diff(start_time, stop_time);

    MPI_Barrier(MPI_COMM_WORLD);
    dfs_level_cutoff++;
  }
  if(store_embeddings && store_in_file)
    dbio::remove_vmap_file(vmap_filename_prefix.c_str());
  //dbio::remove_vmap_file(vmap_filename_prefix);

  if(rank == 0)
    gs_handler.print_profile_info();

  for(int i = 0; i < 5; i++) {
    double sum = 0.0;
    MPI_Reduce(&profile_timings[profiles[i]], &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //cout<<"rank "<<rank<<" profile = "<<profiles[i]<<" time = "<<profile_timings[profiles[i]]<<endl;
    if(rank == 0)
      INFO(*(graph_miner::logger), "AVG time for " << profiles[i] << " = " << (sum / numtasks));
  }
}

void graph_miner_mpi_omp_part::compute_exact_support(){

  timeval start_time, stop_time;
  double profile_timing = 0.0;

  //gather all globally frequent dfs codes
  std::vector<types::DFSCode> globally_frequent;
  gs_handler.get_all_globally_estimated_patterns(globally_frequent);

  DEBUG(*(graph_miner::logger), "The globally estimated patterns for false positives test= " << globally_frequent.size());

  //all ranks have the same order of the frequent patterns
  std::sort(globally_frequent.begin(), globally_frequent.end());

  //for 4 ranks and 2 partitions , groups (0,1) (2,3)
  MPI_Comm rank_group;
  int group_id = rank / num_partitions;
  MPI_Comm_split(MPI_COMM_WORLD, group_id, rank, &rank_group);

  int num_groups = numtasks / num_partitions;
  int patterns_per_group =  (int) ceil( globally_frequent.size() * 1.0 / num_groups );
  int start_index = group_id * patterns_per_group;
  int end_index = start_index + patterns_per_group - 1;
  if(end_index > globally_frequent.size() - 1)
    end_index = globally_frequent.size() - 1;

  DEBUG(*(graph_miner::logger), "rank " << rank << " group id " << group_id << " total frequent " << globally_frequent.size() << " start_index " << start_index << " end_index " << end_index);

  const int MAX_PATTERNS = 200; //send as many patterns as possible, too many may decrease performance

  //marks whether start_index + i pattern is a false positive or not, 1 for yes, 0 for no , also support value
  int *false_positive_marker = 0;
  int false_positive_marker_size = 2;

  if(start_index <= end_index) {
    false_positive_marker_size = 2 + (end_index - start_index + 1);
    false_positive_marker = new int[false_positive_marker_size];
    false_positive_marker[0] = start_index;
    false_positive_marker[1] = end_index;

    //handle patterns in parts
    int k = start_index;

    while(k <= end_index) {

      std::map<string, std::vector<std::vector<std::set<int> > > > dfscode_vertex_mappings;
      int num_patterns = 0, total_vertex_map_size = 0;

      if ( k + MAX_PATTERNS > end_index) {
        num_patterns = end_index - k + 1;
      }else{
        num_patterns = MAX_PATTERNS;
      }

      int *vertex_map_sizes = new int[num_patterns];

      gettimeofday(&start_time, 0);

      int numthreads = threads_per_rank;
      if(store_embeddings)
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
          dfscode = globally_frequent.at(k + i);
          //}catch(std::out_of_range o){
          //	 std::cout<<o.what()<<std::endl;
          //}
          string dfs_str = dfscode.to_string();

          int vmap_size = 1;
          for(int m = 0; m <dfscode.size(); m++) {
            if(dfscode[m].is_forward())
              vmap_size++;
          }

          std::vector<std::set<int> > vertex_maps(vmap_size);
          //std::vector<std::set<int> > vertex_maps;
          bool found = false;
          if(store_embeddings) {
            found = gs_handler.get_stored_vertex_mappings(dfscode, store_in_file, vertex_maps);
          }
          if(found == false)
            get_vertex_maps(dfscode, vertex_maps);

          std::vector<std::vector<std::set<int> > > global_ids_by_partitions(num_partitions, std::vector<std::set<int> >(vertex_maps.size() ));
          omp_set_lock(&lock);
          dfscode_vertex_mappings[dfs_str] = global_ids_by_partitions;
          omp_unset_lock(&lock);
          vertex_map_sizes[i] = vertex_maps.size();
          //total_vertex_map_size += vertex_maps.size();
          try{
            for(int j = 0; j < vertex_maps.size(); j++) {
              for( std::set<int>::iterator it = vertex_maps.at(j).begin(); it != vertex_maps.at(j).end(); it++) {
                int local_vid = *it;
                //omp_set_lock(&lock);
                //dfscode_vertex_mappings[dfs_str][graph[local_vid].orig_part_id][j].insert(graph[local_vid].global_vid);
                //dfscode_vertex_mappings.at(dfs_str).at(graph[local_vid].orig_part_id).at(j).insert(graph[local_vid].global_vid);
                dfscode_vertex_mappings[dfs_str].at(graph[local_vid].orig_part_id).at(j).insert(graph[local_vid].global_vid);
                //omp_unset_lock(&lock);
                //INFO(*(graph_miner::logger), " rank "<<rank<<" local vid "<<local_vid<<" owner rank = "<<graph[local_vid].orig_part_id);
              }
            }
          } catch(std::out_of_range o) {
            std::cout << o.what() << std::endl;
          }
        }                         //for
      }                   // omp parallel
                          //}

      gettimeofday(&stop_time, 0);
      profile_timing += utils::get_time_diff(start_time, stop_time);

      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DEBUG(*(graph_miner::logger), "Processing number of patterns " << num_patterns);

      for(int j = 0; j < num_patterns; j++)
        total_vertex_map_size += vertex_map_sizes[j];
      delete [] vertex_map_sizes;

      int *send_lengths = new int[num_partitions];           //number of neighbor requests per rank
      int *recv_lengths = new int[num_partitions];
      int total_send_length = 0, total_recv_length = 0;

      int *send_buf = 0;
      int *recv_buf = 0;
      int *cum_recv_lengths = new int[num_partitions];
      int offset = 0;

      //#pragma omp parallel for
      for(int j = 0; j < num_partitions; j++) {
        send_lengths[j] = 0;
        if(j != rank % num_partitions) {
          for(int i = 0; i < num_patterns; i++) {
            types::DFSCode dfscode = globally_frequent[k + i];
            string dfs_str = dfscode.to_string();

            send_lengths[j] += dfscode_vertex_mappings[dfs_str][j].size();                                             // number of mappings for pattern vertices
            for(int l = 0; l < dfscode_vertex_mappings[dfs_str][j].size(); l++)
              send_lengths[j] += dfscode_vertex_mappings[dfs_str][j][l].size();
          }
          total_send_length += send_lengths[j];
        }
      }

      MPI_Alltoall(send_lengths, 1, MPI_INT, recv_lengths, 1, MPI_INT, rank_group);
      DEBUG(*(graph_miner::logger)," send lengths " << utils::print_array(send_lengths, num_partitions) << " recv lengths " << utils::print_array(recv_lengths, num_partitions));

      //#pragma omp parallel for
      for(int j = 0; j < num_partitions; j++) {
        total_recv_length += recv_lengths[j];
      }

      //Now allocate send_buf and recv_buf;
      send_buf = new int[total_send_length];
      recv_buf = new int[total_recv_length];

      offset = 0;
      for(int j = 0; j < num_partitions; j++) {
        if(j != rank % num_partitions) {
          for(int i = 0; i < num_patterns; i++) {
            types::DFSCode dfscode = globally_frequent[k + i];
            string dfs_str = dfscode.to_string();
            //stringstream ss;
            //ss <<"DFS CODE "<<dfs_str<<" partition "<<j<<" external vertex mappings ";
            for(int l = 0; l < dfscode_vertex_mappings[dfs_str][j].size(); l++) {
              //ss<< "vertex id = "<<l<<" : ";
              send_buf[offset++] = dfscode_vertex_mappings[dfs_str][j][l].size();
              for(std::set<int>::iterator it  = dfscode_vertex_mappings[dfs_str][j][l].begin(); it!= dfscode_vertex_mappings[dfs_str][j][l].end(); it++) {
                send_buf[offset++] =  *it;
                //ss<< *it<< " ";
              }
            }
            //INFO(*(graph_miner::logger), ss.str());
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

      DEBUG(*(graph_miner::logger), " send_buf " << utils::print_array(send_buf, total_send_length) << " recv_buf " << utils::print_array(recv_buf, total_recv_length));

      //Now process received data, compute actual local support
      offset = 0;
      for(int j = 0; j < num_partitions; j++) {
        if(j != rank % num_partitions) {
          for(int i = 0; i < num_patterns; i++) {
            types::DFSCode dfscode = globally_frequent[k + i];
            string dfs_str = dfscode.to_string();

            for(int l = 0; l < dfscode_vertex_mappings[dfs_str][rank % num_partitions].size(); l++) {
              int len = recv_buf[offset++];
              for(int m = 0; m < len; m++) {
                dfscode_vertex_mappings[dfs_str][rank % num_partitions][l].insert(recv_buf[offset++]);
              }
            }
          }
        }
      }

      //Allreduce
      int *sendcounts = new int[total_vertex_map_size];
      int *recvcounts = new int[total_vertex_map_size];

      offset = 0;
      for(int i = 0; i < num_patterns; i++) {
        types::DFSCode dfscode = globally_frequent[k + i];
        string dfs_str = dfscode.to_string();
        for(int l = 0; l < dfscode_vertex_mappings[dfs_str][rank % num_partitions].size(); l++) {                        //vertex map size
          sendcounts[offset++] = dfscode_vertex_mappings[dfs_str][rank % num_partitions][l].size();
        }
      }

      MPI_Allreduce(sendcounts, recvcounts, total_vertex_map_size, MPI_INT, MPI_SUM, rank_group);

      //compute support
      offset = 0;

      for(int i = 0; i < num_patterns; i++) {
        types::DFSCode dfscode = globally_frequent[k + i];
        string dfs_str = dfscode.to_string();
        int min = 0x7FFFFFFF;
        for(int l = 0; l < dfscode_vertex_mappings[dfs_str][rank % num_partitions].size(); l++) {                        //vertex map size
          if(min > recvcounts[offset])
            min = recvcounts[offset];
          offset++;
        }


        if(min < minimal_support) {
          if(rank == 0)
            num_false_positives++;
          DEBUG(*(graph_miner::logger), "DFS code" << dfscode.to_string() << " has support " << min << " is not globally frequent");
          //Remove if present in local list
          gs_handler.remove_globally_frequent_dfscode(dfscode);
          false_positive_marker[2 + k - start_index + i] = 1;
        }else{
          //update current support
          gs_handler.set_dfscode_support(dfscode, min);
          false_positive_marker[2 + k - start_index + i] = 0;
          if(rank == 0) {
            //report(dfscode, min);
            DEBUG(*(graph_miner::logger), "Globally frequent = " << dfscode.to_string() << " : support = " << min);
          }
        }
      }

      delete [] send_lengths;
      delete [] recv_lengths;

      if(send_buf)
        delete [] send_buf;
      if(recv_buf)
        delete [] recv_buf;

      delete [] sendcounts;
      delete [] recvcounts;
      //delete [] send_buf2;
      //delete [] recv_buf2;
      delete [] cum_recv_lengths;

      k += MAX_PATTERNS;
    }

  }else{
    false_positive_marker = new int[2];
    false_positive_marker[0] = start_index;
    false_positive_marker[1] = end_index;
    false_positive_marker_size = 2;
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

    DEBUG(*(graph_miner::logger), " Lengths received: " << utils::print_array(lens, num_mirrors));

    int *displs = new int[num_mirrors];
    displs[0] = 0;
    for(int i = 1; i< num_mirrors; i++)
      displs[i] = displs[i - 1] + lens[i - 1];

    MPI_Allgatherv(false_positive_marker, false_positive_marker_size, MPI_INT, recv_buf2, lens, displs, MPI_INT, rank_mirror);
    delete [] displs;

    //INFO(*(graph_miner::logger), " Sent "<<utils::print_array(false_positive_marker, false_positive_marker_size) );
    //INFO(*(graph_miner::logger), " Received "<<utils::print_array(recv_buf2, total_recv_length));

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
            types::DFSCode dfscode = globally_frequent[j];
            int is_false_pos = recv_buf2[offset++];
            //int sup = recv_buf2[offset++];
            DEBUG(*(graph_miner::logger), "DFS Code " << globally_frequent[j].to_string() << " is false pos" << is_false_pos << " sup " << sup);
            if(is_false_pos) {
              if(rank == 0)
                num_false_positives++;
              //Remove if present in local list
              gs_handler.remove_globally_frequent_dfscode(dfscode);
            }else{
              //update current support
              //gs_handler.set_dfscode_support(dfscode, sup);
              DEBUG(*(graph_miner::logger), "DFS code" << dfscode.to_string() << " has support " << min << " is not globally frequent");
            }
          }
        }

      }
    }

    delete [] lens;
    delete [] recv_buf2;
  }

  delete [] false_positive_marker;

  //int cfreq_rank = gs_handler.get_globally_frequent(), cumulative_freq;
  //MPI_Reduce(&cfreq_rank, &cumulative_freq, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank == 0) {
    //INFO(*(graph_miner::logger), "FALSE POSITIVES = " << num_false_positives<< " cumulative frequent patterns ="<<cumulative_freq);
    INFO(*(graph_miner::logger), "FALSE POSITIVES = " << num_false_positives);
  }
  INFO(*(graph_miner::logger), "rank = " << rank << " dfs level = " << dfs_level_cutoff << " regeneration time = " << profile_timing);
}


int graph_miner_mpi_omp_part::get_num_neighbors(int src, int global_vid){

  int num_neighbors = 0;
  if(graph.global_local_id_map.count(global_vid) > 0) {
    int local_vid = graph.global_local_id_map[global_vid];
    for(int i = 0; i < graph[local_vid].edge.size(); i++)
      if(graph[graph[local_vid].edge[i].to].orig_part_id != src)
        num_neighbors++;
  }

  return num_neighbors;
}

int graph_miner_mpi_omp_part::get_num_neighbors(int src, int global_vid, std::set<int> &exclusions){

  int num_neighbors = 0;
  if(graph.global_local_id_map.count(global_vid) > 0) {
    int local_vid = graph.global_local_id_map[global_vid];
    for(int i = 0; i < graph[local_vid].edge.size(); i++)
      if(graph[graph[local_vid].edge[i].to].orig_part_id != src
         && exclusions.find(graph[graph[local_vid].edge[i].to].global_vid) == exclusions.end() )
        num_neighbors++;
  }

  return num_neighbors;
}

int graph_miner_mpi_omp_part::prepare_neighbor_data_to_send(int src, int global_vid, int *buffer, int size){
  int pos = 0;
  if(graph.global_local_id_map.count(global_vid) > 0) {
    pos = graph.serialize_neighbors_for_partition(src, global_vid, buffer, size);
    DEBUG(*(graph_miner::logger), "sending neighbor data size = " << pos);
  }
  return pos;
}

void graph_miner_mpi_omp_part::process_received_neighbor_data(int *buffer, int size){
  int global_vid;
  graph.deserialize_neighbors_for_partition(rank, buffer, size, global_vid);
}

void graph_miner_mpi_omp_part::process_received_neighbor_data(int *buffer, int size, int num_neighbors){
  graph.deserialize_multiple_neighbors_for_partition(rank % num_partitions, buffer, size, num_neighbors);
}

/*
 * Invoked when purge_external_neighbors flag is turned on
 * Tries to clean up the non used pseudo locals
 */
void graph_miner_mpi_omp_part::remove_non_used_pseudo_locals(){

  int LRU_levels_to_consider = 2;
  if(dfs_level_cutoff - LRU_levels_to_consider  <= 0)
    return;

  std::set<int> used_pseudo_locals;       //local ids

  //store in set local_ids of all used pseudo locals for the last "LRU_levels_to_consider" iterations
  for(int l = dfs_level_cutoff; l>= dfs_level_cutoff - LRU_levels_to_consider; l--) {
    for(set<int>::iterator it = pseudo_locals_on_extensions_per_level[l].begin();
        it != pseudo_locals_on_extensions_per_level[l].end(); it++) {
      // it is global id
      int local_id = graph.get_local_vid(*it);
      used_pseudo_locals.insert(*it);
    }
  }

  //sync pseudo locals with mirror ranks
  if(numtasks > num_partitions)
    union_used_pseudo_locals_in_rank_group(used_pseudo_locals);

  std::vector<int> to_delete_vertices;        //local ids
  std::set<int> to_downgrade_vertices;       //global ids

  while(true) {
    for(int i = graph.max_local_vid + 1; i< graph.size(); i++) {
      //if external vertex
      if(graph.is_external(i)) {
        bool only_connected_to_not_used_pl = true;
        std::set<int> neighbor_global_ids;                         //global ids of neighbors
        //iterate through all it's neighbors
        for(int k = 0; k < graph[i].edge.size(); k++) {
          int to_vertex = graph[i].edge[k].to;
          neighbor_global_ids.insert(graph[to_vertex].global_vid);
          //if the neighbor is not a "not used" pseudo local, i cannot be deleted
          if(!graph.is_pseudo_local(to_vertex) || used_pseudo_locals.count(graph[to_vertex].global_vid) > 0) {
            only_connected_to_not_used_pl = false;
            break;
          }
        }

        if(only_connected_to_not_used_pl) {
          to_delete_vertices.push_back(i);
          for(std::set<int>::iterator it = neighbor_global_ids.begin(); it != neighbor_global_ids.end(); ++it)
            to_downgrade_vertices.insert(*it);
        }

      }
    }

    //if delete list empty break
    if(to_delete_vertices.size()  == 0)
      break;

    //Now downgrade the pseudo locals
    for(std::set<int>::iterator it = to_downgrade_vertices.begin(); it != to_downgrade_vertices.end(); ++it) {
      int local_id = graph.get_local_vid(*it);
      graph[local_id].vertex_part_id = graph[local_id].orig_part_id;
    }

    //And delete the not needed externals from the graph
    if(to_delete_vertices.size() > 0)
      graph.delete_vertices(to_delete_vertices);

    //now the downgraded vertices are externals
    //iterate through these externals and remove inconsistency ( edge between external - external)
    for(std::set<int>::iterator it = to_downgrade_vertices.begin(); it != to_downgrade_vertices.end(); it++) {
      //*it is global id
      int local_vid = graph.get_local_vid(*it);
      if(local_vid == -1) continue;
      //iterate through all it's neighbors
      for(int k = 0; k < graph[local_vid].edge.size(); k++) {
        int to_vertex = graph[local_vid].edge[k].to;
        //remove inconsistency in graph
        //if the neighbor is an external vertex, delete edge between two externals
        if(graph.is_external(to_vertex)) {
          graph.delete_edge(local_vid, to_vertex);
        }
      }
    }

    INFO(*(graph_miner::logger), " Pseudo locals on extensions for last two iterations = " << used_pseudo_locals.size()
         << " deleted external " << to_delete_vertices.size() << " downgraded " << to_downgrade_vertices.size() );

    to_downgrade_vertices.clear();
    to_delete_vertices.clear();
  }      //while

  //clean up the unused level for the next iteration
  if(dfs_level_cutoff - 2 > 0)
    pseudo_locals_on_extensions_per_level[dfs_level_cutoff - 2].clear();
}

} // namespace GRAPH_MINER
