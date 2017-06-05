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

#include <graph_miner_mpi_omp_hybrid.hpp>
#include <utils.hpp>
#include <iterator>
#include <stdlib.h>
#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <math.h>
#include <omp.h>

using namespace std;

namespace GRAPH_MINER {

graph_miner_mpi_omp_hybrid::graph_miner_mpi_omp_hybrid(int th_per_rank) : graph_miner()
{
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  task_split_threshold = 2;
  is_working = true;
  donor_thread = 0;

  DEBUG(*(graph_miner::logger), " omp_max_threads = " << omp_get_max_threads());
  if(th_per_rank < omp_get_max_threads()) {
    threads_per_rank = omp_get_max_threads();
  }else{
    threads_per_rank = th_per_rank;
  }

  for(int i = 0; i < numtasks; i++) {
    std::vector<std::vector<std::deque<types::DFS> > > tmp;
    global_shared_queue.push_back(tmp);
    for(int j = 0; j < threads_per_rank; j++) {
      std::vector<std::deque<types::DFS> > tmp2;
      global_shared_queue[i].push_back(tmp2);
    }
  }

  for(int i = 0; i<threads_per_rank; i++) {

    frequent_patterns_count.push_back(0);
    task_split_level.push_back(0);     // splits are done from the top of the tree first
    embeddings_regeneration_level.push_back(0);
    current_dfs_level.push_back(0);
    thread_is_working.push_back(false);
    types::DFSCode dfscode;
    DFS_CODE_V.push_back(dfscode);
    DFS_CODE_IS_MIN_V.push_back(dfscode);
    types::Graph gr;
    GRAPH_IS_MIN_V.push_back(gr);
    std::vector<std::deque<types::DFS> > tmp;
    dfs_task_queue.push_back(tmp);
    dfs_task_queue_shared.push_back(tmp);
    omp_lock_t l;
    omp_init_lock(&l);
    qlock.push_back(&l);

  }

  //MPI_Comm comm_object;
  //int membershipKey = rank / numtasks;
  //MPI_Comm_split(MPI_COMM_WORLD, membershipKey, rank, &comm_object);
  //init_dynamic_load_balancing(rank, numtasks, comm_object);
  init_dynamic_load_balancing(rank, numtasks, MPI_COMM_WORLD);

  init_dynamic_threads_load_balancing(rank, threads_per_rank);

}


graph_miner_mpi_omp_hybrid::~graph_miner_mpi_omp_hybrid(){

  for(int i = 0; i<threads_per_rank; i++) {
    omp_destroy_lock(qlock[i]);
  }

  qlock.empty();
}


void graph_miner_mpi_omp_hybrid::set_task_split_threshold(int thr){
  task_split_threshold = thr;
}

void graph_miner_mpi_omp_hybrid::set_load_balance_interval(int i){
  load_balance_interval = i;
}

void graph_miner_mpi_omp_hybrid::set_num_threads(int t){
  omp_set_num_threads(t);
}


int graph_miner_mpi_omp_hybrid::get_frequent_patterns_count(){
  int total = 0;
  for(int i = 0; i<threads_per_rank; i++)
    total += frequent_patterns_count[i];
  return total;
}

void graph_miner_mpi_omp_hybrid::print_threads_frequent_patterns(){
  stringstream ss;
  for(int i = 0; i < threads_per_rank; i++)
    ss << " tid = " << i << " patterns = " << frequent_patterns_count[i];
  INFO(*(graph_miner::logger), ss.str());
}

void graph_miner_mpi_omp_hybrid::print_threads_status(int thread_id){
  stringstream ss;
  ss << " From thread = " << thread_id;
  for(int i = 0; i < threads_per_rank; i++) {
    if(thread_is_working[i] == true)
      ss << " tid = " << i << " active ";
    else
      ss << " tid = " << i << " idle ";
  }
  INFO(*(graph_miner::logger), ss.str());
}

//support function for a single large graph, computes the minimum count of a node in the embeddings
unsigned int
graph_miner_mpi_omp_hybrid::support(Projected &projected, Thread_private_data &gprv)
{
  std::map<unsigned int, map<unsigned int, unsigned int> > node_id_counts;
  //Print DFS code
  //for(int i = 0; i < DFS_CODE.size(); i++)
  //  std::cout << DFS_CODE[i].to_string();
  //std::cout << endl;

  int thread_id = gprv.thread_id; //omp_get_thread_num();
  //iterated through the all the embeddings
  for(Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {

    PDFS *em = &(*cur);
    int dfsindex = gprv.DFS_CODE.size() - 1;
    while(em != NULL) {
      //std::cout<<DFS_CODE[dfsindex].to_string() <<endl;
      //print embedding
      //std::cout<<em->to_string()<<endl;
      if(gprv.DFS_CODE[dfsindex].to > gprv.DFS_CODE[dfsindex].from) {    //forward edge
        node_id_counts[gprv.DFS_CODE[dfsindex].to][em->edge->to]++;
      }
      if(!em->prev) {
        node_id_counts[gprv.DFS_CODE[dfsindex].from][em->edge->from]++;
      }
      em = em->prev;
      dfsindex--;
    }
    //std::cout<<endl;
  }

  unsigned int min = 0xffffffff;
  for(std::map<unsigned int, map<unsigned int, unsigned int> >::iterator it = node_id_counts.begin(); it != node_id_counts.end(); it++) {
    if((it->second).size() < min)
      min = (it->second).size();
  }
  if(min == 0xffffffff) {
    //std::cout << gprv.DFS_CODE.to_string()<<" support: " << min << endl;
    min = 0;
  }
  return min;

}

/* Recursive subgraph mining function (similar to subprocedure 1
 * Subgraph_Mining in [Yan2002]).
 */
void graph_miner_mpi_omp_hybrid::project(Projected &projected, int dfs_level, Thread_private_data &gprv)
{

  // Check if the pattern is frequent enough.
  unsigned int sup = support(projected, gprv);
  if(sup < minimal_support) return;

  int thread_id = gprv.thread_id; //omp_get_thread_num();

  // The minimal DFS code check is more expensive than the support check,
  // hence it is done now, after checking the support.
  if(is_min(gprv) == false) {
    return;
  } else {

  }

  DEBUG(*logger, "DFS level = " << dfs_level);
  //omp_set_lock(&lock);
  //INFO(*(graph_miner::logger), "thread " << thread_id <<" executing project for code: " << gprv.DFS_CODE.to_string() << "; support: " << sup);
  //omp_unset_lock(&lock);

  // Output the frequent substructure
  report(projected, sup);
  frequent_patterns_count[thread_id]++;

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
  int maxtoc = gprv.DFS_CODE[rmpath[0]].to;

  Projected_map3 new_fwd_root;
  Projected_map2 new_bck_root;
  types::EdgeList edges;

  gprv.current_dfs_level = dfs_level;

  // Enumerate all possible one edge extensions of the current substructure.
  for(unsigned int n = 0; n < projected.size(); ++n) {

    unsigned int id = projected[n].id;
    PDFS *cur = &projected[n];
    History history(graph, cur);

    // XXX: do we have to change something here for directed edges?

    // backward
    for(int i = (int)rmpath.size() - 1; i >= 1; --i) {
      Edge *e = get_backward(graph, history[rmpath[i]], history[rmpath[0]], history);
      if(e) {
        new_bck_root[gprv.DFS_CODE[rmpath[i]].from][e->elabel].push(id, e, cur);
      }
    }

    // pure forward
    // FIXME: here we pass a too large e->to (== history[rmpath[0]]->to
    // into get_forward_pure, such that the assertion fails.
    //
    // The problem is:
    // history[rmpath[0]]->to > graph_database[id].size()
    if(get_forward_pure(graph, history[rmpath[0]], minlabel, history, edges)) {
      for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it) {
        new_fwd_root[maxtoc][(*it)->elabel][graph[(*it)->to].label].push(id, *it, cur);
      }
    }

    // backtracked forward
    for(int i = 0; i < (int)rmpath.size(); ++i) {
      if(get_forward_rmpath(graph, history[rmpath[i]], minlabel, history, edges)) {
        for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it) {
          new_fwd_root[gprv.DFS_CODE[rmpath[i]].from][(*it)->elabel][graph[(*it)->to].label].push(id, *it, cur);
        } // for it
      } // if
    } // for i
  } // for n

  std::deque<types::DFS> tmp;

  if(gprv.dfs_task_queue.size() <= dfs_level) {
    gprv.dfs_task_queue.push_back(tmp);
  }

  // Test all extended substructures.
  // backward
  for(Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
    for(Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {

      DFS dfs(maxtoc, to->first, -1, elabel->first, -1);
      gprv.dfs_task_queue[dfs_level].push_back(dfs);

    }
  }

  // forward
  for(Projected_riterator3 from = new_fwd_root.rbegin();
      from != new_fwd_root.rend(); ++from) {
    for(Projected_iterator2 elabel = from->second.begin();
        elabel != from->second.end(); ++elabel) {
      for(Projected_iterator1 tolabel = elabel->second.begin();
          tolabel != elabel->second.end(); ++tolabel) {

        DFS dfs(from->first, maxtoc + 1, -1, elabel->first, tolabel->first);
        gprv.dfs_task_queue[dfs_level].push_back(dfs);

      }
    }
  }

  while(gprv.dfs_task_queue[dfs_level].size() > 0) {

    if(numtasks > 1 && thread_id == 0)
      load_balance(gprv);

    if(num_threads > 1)
      threads_load_balance(gprv);


    DFS dfs = gprv.dfs_task_queue[dfs_level].front();
    gprv.dfs_task_queue[dfs_level].pop_front();

    DEBUG(*(graph_miner::logger), "popped dfs = " << dfs.to_string() );

    gprv.current_dfs_level = dfs_level;

    gprv.DFS_CODE.push(dfs.from, dfs.to, dfs.fromlabel, dfs.elabel, dfs.tolabel);

    if(dfs.is_backward())
      project(new_bck_root[dfs.to][dfs.elabel], dfs_level + 1, gprv);      //Projected (PDFS vector): each entry contains graph id 0, edge pointer, null PDFS
    else
      project(new_fwd_root[dfs.from][dfs.elabel][dfs.tolabel], dfs_level + 1, gprv);      //Projected (PDFS vector): each entry contains graph id 0, edge pointer, null PDFS

    gprv.DFS_CODE.pop();

  }

  return;
}

/*
 * Regnerate the embeddings for the dfs codes obtained from other processors
 */

void graph_miner_mpi_omp_hybrid::regenerate_embeddings(Projected &projected, int dfs_level, Thread_private_data &gprv)
{
  // We don't need to check if the pattern is frequent or minimal

  DEBUG(*(graph_miner::logger), "DFS level inside regenerate embeddings = " << dfs_level << " queue size = " << dfs_task_queue[dfs_level].size());

  //not necessary though, as task split is not done while regenerating embeddings
  //current_dfs_level = dfs_level + 1;

  int thread_id = gprv.thread_id; //omp_get_thread_num();

  //iterate for all in the task_queue
  for(int i = 0; gprv.dfs_task_queue[dfs_level].size() > 0; i++) {

    //std::cout << "dfs level = " << dfs_level <<std::endl;
    gprv.current_dfs_level = dfs_level;
    if(thread_id == 0)
      load_balance(gprv);
    if(num_threads > 1)
      threads_load_balance(gprv);

    types::DFS dfs = gprv.dfs_task_queue[dfs_level].front();
    gprv.dfs_task_queue[dfs_level].pop_front();

    gprv.DFS_CODE.push(dfs.from, dfs.to, dfs.fromlabel, dfs.elabel, dfs.tolabel);

    //omp_set_lock(&lock);
    //INFO(*(graph_miner::logger), "thread " << thread_id << " *****regenerating embeddings for code: " << gprv.DFS_CODE.to_string() );
    //omp_unset_lock(&lock);

    Projected new_root;

    for(unsigned int n = 0; n < projected.size(); ++n) {

      unsigned int id = projected[n].id;
      PDFS *cur = &projected[n];
      History history(graph, cur);

      if(dfs.is_backward() ) {
        Edge *e = get_backward(graph, gprv.DFS_CODE, history);
        if(e)
          new_root.push(id, e, cur);
      }else{
        types::EdgeList edges;
        if(get_forward(graph, gprv.DFS_CODE, history, edges)) {
          for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it) {
            new_root.push(id, *it, cur);
          }
        }
      }
    }

    if( gprv.embeddings_regeneration_level > dfs_level ) {
      regenerate_embeddings(new_root, dfs_level + 1, gprv);
    }else{
      //regeneration of embeddings ended
      //now perform regular extensions with project function
      //reset embeddings_regeneration_level
      //embeddings_regeneration_level = 0;
      project(new_root, dfs_level + 1, gprv);
    }

    gprv.DFS_CODE.pop();
  }

  return;

}

void graph_miner_mpi_omp_hybrid::run_intern(void)
{

  types::EdgeList edges;
  Projected_map3 root;
  int single_edge_dfscodes = 0;

  for(unsigned int from = 0; from < graph.size(); ++from) {
    if(get_forward_root(graph, graph[from], edges)) {   // get the edge list of the node g[from] in graph g
      for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it) {
        //embeddings with a single edge
        if(root.count(graph[from].label) == 0 || root[graph[from].label].count((*it)->elabel) == 0 || root[graph[from].label][(*it)->elabel].count(graph[(*it)->to].label) == 0) {
          single_edge_dfscodes++;
          DEBUG(*logger, "single edge DFS code : (0,1," << graph[from].label << "," << (*it)->elabel << "," << graph[(*it)->to].label << ")" );
        }
        root[graph[from].label][(*it)->elabel][graph[(*it)->to].label].push(0, *it, 0);          //projected (PDFS vector) entry: graph id (always 0 for single graph), edge pointer and null PDFS
      }  //for
    }   // if
  }   // for from
  //} // for id

  Thread_private_data gprv;

   #pragma omp parallel num_threads(threads_per_rank) private(gprv)
  {

    int thread_id = omp_get_thread_num();
    gprv.thread_id = thread_id;
    //gprv.load_balance_interval = 5; //load_balance_interval;
    gprv.current_dfs_level = 0;
    gprv.task_split_level = 0;
    gprv.embeddings_regeneration_level = 0;
    //gprv.frequent_patterns_count = 0;

    Projected emb_prv;

    int dfscodes_per_rank =  (int) ceil((single_edge_dfscodes * 1.0) / numtasks);
    int dfscodes_per_thread =  (int) ceil((dfscodes_per_rank * 1.0) / threads_per_rank);

    int start_index = rank * dfscodes_per_rank + thread_id * dfscodes_per_thread;
    int end_index = start_index + dfscodes_per_thread - 1;

    if (end_index >= (rank + 1) * dfscodes_per_rank)
      end_index = (rank + 1) * dfscodes_per_rank - 1;

    if (end_index > single_edge_dfscodes - 1)
      end_index = single_edge_dfscodes - 1;


    if(start_index <= end_index) { //has some work
      omp_set_lock(&lock);
      thread_is_working[thread_id] = true;
      omp_unset_lock(&lock);
    }

    DEBUG(*(graph_miner::logger), "thread id " << thread_id << "start index = " << start_index << " , end index = " << end_index << endl);

    std::deque<types::DFS> tmp;
    gprv.dfs_task_queue.push_back(tmp);

    int index = 0;
    for(Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end(); ++fromlabel) {
      for(Projected_iterator2 elabel = fromlabel->second.begin(); elabel != fromlabel->second.end(); ++elabel) {
        for(Projected_iterator1 tolabel = elabel->second.begin();
            tolabel != elabel->second.end(); ++tolabel) {

          if( index >= start_index && index <= end_index ) {
            // Build the initial two-node graph.  It will be grownrecursively within project.

            DFS dfs(0, 1, fromlabel->first, elabel->first, tolabel->first);
            gprv.dfs_task_queue[0].push_back(dfs);
          }
          index++;

        } // for tolabel
      } // for elabel
    } // for fromlabel

    //TRACE((*graph_miner::logger), "size = " << dfs_task_queue[thread_id][0].size() );

    while(computation_end == false) {

      //TRACE(*(graph_miner::logger),"processor "<<rank<<" thread "<< thread_id <<" is idle = " << thread_is_working[thread_id]);

      if(thread_working(gprv) == false) {

        if(num_threads > 1 && computation_end == false)            //&& thread_is_working[thread_id] == false)
          threads_load_balance(gprv);

        if(thread_id == 0) {

          if(all_threads_idle() == true) {
            is_working = false;
            if(numtasks > 1) {
              bool send_req = false;
              //thread 0 collects all work requests from all threads before making a global work request (rank)
              if(all_idle_work_requests_collected()) {
                send_req = true;
                reset_idle_work_requests_counter();
              }
              load_balance(gprv, send_req);
            }else{
              computation_end = true;                             //single process, computation ended
            }
          }

        }


      }else{
        //this is done in process_received_data, so not required here
        //thread_is_working[thread_id] = true;

        if(numtasks > 1 && thread_id == 0)
          load_balance(gprv);

        if(num_threads > 1)
          threads_load_balance(gprv);

        DFS dfs = gprv.dfs_task_queue[0].front();

        gprv.dfs_task_queue[0].pop_front();

        //DFS dfs = dfs_task_queue[thread_id][0][i];
        //omp_set_lock(&lock);
        DEBUG(*(graph_miner::logger), "thread " << thread_id << " popped dfs = " << dfs.to_string() << " regen level = " << gprv.embeddings_regeneration_level );
        //omp_unset_lock(&lock);
        //omp_set_lock(&lock);
        gprv.DFS_CODE.push(0, 1, dfs.fromlabel, dfs.elabel, dfs.tolabel);
        //omp_unset_lock(&lock);

        gprv.current_dfs_level = 1;
        //emb_prv = root[dfs.fromlabel][dfs.elabel][dfs.tolabel];
        //TRACE(*(graph_miner::logger), "embeddings regeneration level = " << embeddings_regeneration_level);
        if(gprv.embeddings_regeneration_level < 1)
          //project(emb_prv, 1, gprv); //Projected (PDFS vector): each entry contains graph id 0, edge pointer, null PDFS
          project(root[dfs.fromlabel][dfs.elabel][dfs.tolabel], 1, gprv);                  //Projected (PDFS vector): each entry contains graph id 0, edge pointer, null PDFS
        else
          //regenerate_embeddings(emb_prv, 1, gprv);
          regenerate_embeddings(root[dfs.fromlabel][dfs.elabel][dfs.tolabel], 1, gprv);
        gprv.current_dfs_level = 0;
        //omp_set_lock(&lock);
        gprv.DFS_CODE.pop();
        //omp_unset_lock(&lock);

        if(gprv.dfs_task_queue[0].size() == 0) {
          omp_set_lock(&lock);
          thread_is_working[thread_id] = false;
          omp_unset_lock(&lock);
          gprv.embeddings_regeneration_level = 0;
          embeddings_regeneration_level[thread_id] = 0;
          DEBUG(*(graph_miner::logger),"processor " << rank << " thread " << thread_id << " is idle, has token = " << has_token);
          if(has_token == true)
            DEBUG(*(graph_miner::logger),"processor " << rank << " thread " << thread_id << " token color = " << get_token_color(token_color));
        }

      }
    }

    //omp_destroy_lock(&lock);

    TRACE(*(graph_miner::logger), "Done... rank = " << rank << " thread " << thread_id << " computation ended = " << computation_end_shared[thread_id]);
    if(thread_id == 0) {
      print_threads_frequent_patterns();
      print_threads_status(0);
    }

  }  //pragma omp

} // void graph_miner_mpi_omp_hybrid::run_intern(void)

/* Implemented functions of dynamic_load_balancing class */

bool graph_miner_mpi_omp_hybrid::working(){
  return is_working;
}

void graph_miner_mpi_omp_hybrid::initiate_global_split_work(int requester_rank_id, Thread_private_data &gprv){
  int thread_id = gprv.thread_id;       //omp_get_thread_num();
  DEBUG(*(graph_miner::logger), "thread " << thread_id << " inside initiate_global_split_work");
  //if(thread_id == 0)
  //thread 0 puts work in global shared queue
  thread_split_global_work(requester_rank_id, gprv);
  //thread 0 asks other threads to put work in global shared queue
  dynamic_threads_load_balancing::send_global_split_requests(requester_rank_id, gprv);
}

void graph_miner_mpi_omp_hybrid::complete_global_split_work(int requester_rank_id, Thread_private_data &gprv){

  int MAX_SEND_SIZE = 900;
  //determine the number of dfscodes and the total send length
  int total_send_length = threads_per_rank;        // offsets for threads
  int total_dfscodes = 0;
  for(int i = 0; i < threads_per_rank; i++) {
    int split_level = global_shared_queue[requester_rank_id][i].size() - 1;
    int num_dfs = 0;
    if(split_level >= 0) {
      num_dfs = global_shared_queue[requester_rank_id][i][split_level].size();
      total_send_length += 1 + 5 * ( split_level + num_dfs);                  // task split_level + DFS codes
      total_dfscodes += num_dfs;
    }
    DEBUG(*(graph_miner::logger), " For rank " << requester_rank_id << " thread " << i << "split level " << split_level << " num dfs " << num_dfs << " total_length " << total_send_length);
  }

  int buffer[2];
  buffer[0] = RT_DATA;

  DEBUG(*(graph_miner::logger), "rank " << rank << " thread " << gprv.thread_id << " Sending global split response , send length = " << total_send_length);

  if(total_dfscodes == 0) {
    // send "no work" available
    buffer[1] = 0;             //length 0
    send_other_msg(buffer, 2, requester_rank_id);
  } else{

    //if too big, send in parts
    if(total_send_length > MAX_SEND_SIZE)
      buffer[0] = RT_DATA_PART;

    buffer[1] = total_send_length;             //length 0
    send_other_msg(buffer, 2, requester_rank_id);

    DEBUG(*(graph_miner::logger), "allocating send buf of size = " << total_send_length);
    //serialize work form the global work queue
    int *send_buf = new int[total_send_length];
    DEBUG(*(graph_miner::logger), "succesfully allocated send buf of size = " << total_send_length);
    // data format: first put offsets of all threads data
    // for each thread: split level + DFS codes (prefix + final level DFS candidates)
    int offset = threads_per_rank;
    for(int k = 0; k <threads_per_rank; k++) {
      send_buf[k] = offset;
      int split_level = global_shared_queue[requester_rank_id][k].size() - 1;
      DEBUG(*(graph_miner::logger), "Thread id " << k << " split level = " << split_level << " offset " << offset  );

      int num_dfs = 0;
      if (split_level >= 0) {
        num_dfs = global_shared_queue[requester_rank_id][k][split_level].size();
        send_buf[offset++] = split_level;
        DEBUG(*(graph_miner::logger), "Thread id " << k << " sending split level = " << split_level << " num dfs " << num_dfs  );
      }

      for(int i = 0; i < split_level; i++) {
        global_shared_queue[requester_rank_id][k][i][0].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
        offset += 5;
      }

      for(int j = 0; j < num_dfs; j++) {
        global_shared_queue[requester_rank_id][k][split_level][j].serialize( (char*) (send_buf + offset), sizeof(int) * 5 );
        offset += 5;
      }

      DEBUG(*(graph_miner::logger), "Thread id " << k << " split level = " << split_level );
      //clean up the global shared queue
      global_shared_queue[requester_rank_id][k].clear();
      //if(k == 0 && num_dfs > 0)
      //if(requester_rank_id < rank) my_color = BLACK;

    }

    if(offset !=  total_send_length)
      CRITICAL_ERROR(*(graph_miner::logger), "Problem in serializing data: offset = " << offset << " total length " << total_send_length );
    DEBUG(*(graph_miner::logger), "Sending data " << " total send length = " << total_send_length << " buf = " << utils::print_array(send_buf,total_send_length));

    //send data to the requesting rank
    if(total_send_length <= MAX_SEND_SIZE ) {
      send_data_msg(send_buf, total_send_length, requester_rank_id);
    }else{             // send in parts if data is too big
      offset = 0;
      buffer[0] = RT_DATA;
      while(offset < total_send_length) {
        if(offset + MAX_SEND_SIZE <= total_send_length)
          buffer[1] = MAX_SEND_SIZE;
        else
          buffer[1] = total_send_length - offset;
        send_other_msg(buffer, 2, requester_rank_id);
        send_data_msg(send_buf + offset, buffer[1], requester_rank_id);
        offset += buffer[1];
      }

    }

    //handle Dijkstra's token termination
    if(requester_rank_id < rank) my_color = BLACK;

    delete [] send_buf;
  }

}

void graph_miner_mpi_omp_hybrid::process_received_data(int* buffer, int size, Thread_private_data &gprv){

  //receiver thread id
  int thread_id = gprv.thread_id;       //omp_get_thread_num();

  DEBUG(*(graph_miner::logger), " Received buffer of size = " << size);
  DEBUG(*(graph_miner::logger), " Received buffer " << utils::print_array(buffer, size));

  std::vector<int> thread_has_work;

  for(int k = 0; k < threads_per_rank; k++) {

    //if(k > 0 ) break;
    int offset = buffer[k];
    int received_dfs_level = 0;

    int has_work = 0;
    if( k < threads_per_rank - 1) {
      has_work = (buffer[k + 1] - buffer[k]);
    }else{
      has_work = (size - buffer[k]);
    }

    if(has_work > 0) {
      received_dfs_level = buffer[offset++];
      thread_has_work.push_back(1);
      if(k == 0) {
        //thread 0 can start working
        gprv.embeddings_regeneration_level = received_dfs_level;
      }else{
        embeddings_regeneration_level[k] = received_dfs_level;
      }
    }else{
      thread_has_work.push_back(0);
      continue;
    }

    for(int i = 0; i < received_dfs_level; i++) {

      types::DFS dfs;
      dfs.deserialize((char*) (buffer + offset), 5 * sizeof(int));
      if(k == 0) {                  //thread 0 can put directly into its private queue
        if(gprv.dfs_task_queue.size() < (i + 1) ) {
          std::deque<types::DFS> tmp;
          gprv.dfs_task_queue.push_back(tmp);
        }
        gprv.dfs_task_queue[i].push_front(dfs);
      } else {                   //for other threads put in the global shared queue
        if(dfs_task_queue_shared[k].size() < (i + 1) ) {
          std::deque<types::DFS> tmp;
          dfs_task_queue_shared[k].push_back(tmp);
        }
        dfs_task_queue_shared[k][i].push_front(dfs);
      }
      offset += 5;
      DEBUG( *(graph_miner::logger), "level = " << i << " queue size = " << global_task_queue_shared[k][i].size());
    }

    int num_dfs_received = 0;
    if( k < threads_per_rank - 1) {
      num_dfs_received = (buffer[k + 1] - offset) / 5;
    }else{
      num_dfs_received = (size - offset) / 5;
    }

    DEBUG(*(graph_miner::logger), "Received data : thread " << k << " , split level " << received_dfs_level << " num dfs " << num_dfs_received);

    for(int i = 0; i < num_dfs_received; i++) {
      types::DFS dfs;
      dfs.deserialize((char*) (buffer + offset), 5 * sizeof(int));
      if(k == 0) {                  //thread 0 can put directly into its private queue
        if(gprv.dfs_task_queue.size() < (received_dfs_level + 1) ) {
          std::deque<types::DFS> tmp;
          gprv.dfs_task_queue.push_back(tmp);
        }
        gprv.dfs_task_queue[received_dfs_level].push_front(dfs);
      }else{                   //for other threads put in the global shared queue
        if(dfs_task_queue_shared[k].size() < (received_dfs_level + 1) ) {
          std::deque<types::DFS> tmp;
          dfs_task_queue_shared[k].push_back(tmp);
        }
        dfs_task_queue_shared[k][received_dfs_level].push_front(dfs);
      }
      offset += 5;
      DEBUG(*(graph_miner::logger),  " i = " << i << " received: " << dfs_task_queue_shared[k][received_dfs_level].front());
    }
  }

  dynamic_threads_load_balancing::send_work_response_to_all_threads(gprv, thread_has_work);
  if(thread_has_work[0]) {
    thread_start_working();
  }
  //set process to the ACTIVE state, since there is work for at least one thread
  is_working = true;
}

bool graph_miner_mpi_omp_hybrid::thread_working(Thread_private_data &gprv){

  int thread_id = gprv.thread_id;       //omp_get_thread_num();
  bool th_is_working;
  omp_set_lock(&lock);
  th_is_working = thread_is_working[thread_id];
  omp_unset_lock(&lock);

  return th_is_working;
}

bool graph_miner_mpi_omp_hybrid::all_threads_idle(){

  bool all_idle = true;
  omp_set_lock(&lock);
  for(int i = 0; i<threads_per_rank; i++)
    if(thread_is_working[i] == true) {
      all_idle = false;
      break;
    }
  omp_unset_lock(&lock);
  return all_idle;
}

bool graph_miner_mpi_omp_hybrid::can_thread_split_work(Thread_private_data &gprv){

  int thread_id = gprv.thread_id;      // omp_get_thread_num();

  //return false;
  if(!thread_working(gprv))
    return false;

  //prevents giving work until regeneration is complete
  //if(embeddings_regeneration_level > 0)
  //	return false;

  gprv.task_split_level = 0;

  while(gprv.task_split_level < gprv.current_dfs_level && gprv.dfs_task_queue[gprv.task_split_level].size() < task_split_threshold ) {
    gprv.task_split_level++;
    DEBUG(*(graph_miner::logger), "thread" << thread_id << "task split level changed = " << gprv.task_split_level << " queue size = " << gprv.dfs_task_queue[gprv.task_split_level].size());
  }
  DEBUG(*(graph_miner::logger), "thread" << thread_id << " final task split level = " << gprv.task_split_level << " queue size = " << gprv.dfs_task_queue[gprv.task_split_level].size());

  if( gprv.dfs_task_queue.size() > gprv.task_split_level && gprv.dfs_task_queue[gprv.task_split_level].size() >= task_split_threshold )
    return true;

  return false;

}

void graph_miner_mpi_omp_hybrid::thread_split_work(int requesting_thread, int &length, Thread_private_data &gprv){

  int thread_id = gprv.thread_id;       //omp_get_thread_num();

  DEBUG(*(graph_miner::logger), "inside thread split work : " << gprv.DFS_CODE.to_string() << "  split level " << gprv.task_split_level << ", queue size of level 0 = " << gprv.dfs_task_queue[0].size() );

  for(int i = 0; i < gprv.task_split_level; i++) {
    if(dfs_task_queue_shared[requesting_thread].size() < (i + 1) ) {
      std::deque<types::DFS> tmp;
      //omp_set_lock(&lock);
      dfs_task_queue_shared[requesting_thread].push_back(tmp);
      //omp_unset_lock(&lock);
    }
    //omp_set_lock(&lock);
    dfs_task_queue_shared[requesting_thread][i].push_back(gprv.DFS_CODE[i]);
    //omp_unset_lock(&lock);
  }

  if(dfs_task_queue_shared[requesting_thread].size() < ( gprv.task_split_level + 1) ) {
    std::deque<types::DFS> tmp;
    //omp_set_lock(&lock);
    dfs_task_queue_shared[requesting_thread].push_back(tmp);
    //omp_unset_lock(&lock);
  }

  int num_dfs = gprv.dfs_task_queue[gprv.task_split_level].size() / 2;
  DEBUG(*(graph_miner::logger), "thread " << thread_id << " is splitting work. " << " queue size = " << gprv.dfs_task_queue[gprv.task_split_level].size() << " num dfs to donate = " << num_dfs );

  for(int j = 0; j < num_dfs; j++) {
    //omp_set_lock(&lock);
    types::DFS dfs = gprv.dfs_task_queue[gprv.task_split_level].back();
    dfs_task_queue_shared[requesting_thread][gprv.task_split_level].push_front(dfs);
    DEBUG(*(graph_miner::logger), "thread id " << thread_id << " inside thread split work : " << dfs.to_string() << "  split level " << gprv.task_split_level << ", queue size of level 0 = " << gprv.dfs_task_queue[0].size() );
    gprv.dfs_task_queue[gprv.task_split_level].pop_back();
    //omp_unset_lock(&lock);
  }

  DEBUG(*(graph_miner::logger), "thread " << requesting_thread << " is getting work " << " donor q size = " << gprv.dfs_task_queue[gprv.task_split_level].size() << " requester q size = " << dfs_task_queue[requesting_thread][gprv.task_split_level].size() );
  embeddings_regeneration_level[requesting_thread] = gprv.task_split_level;

  //thread_is_working[requesting_thread] = true;
  length = num_dfs;

}


void graph_miner_mpi_omp_hybrid::thread_process_received_data(Thread_private_data &gprv){

  int thread_id = gprv.thread_id;       //omp_get_thread_num();
  gprv.embeddings_regeneration_level = embeddings_regeneration_level[thread_id];

  int num_dfs = dfs_task_queue_shared[thread_id][gprv.embeddings_regeneration_level].size();
  DEBUG(*(graph_miner::logger), "thread " << thread_id << " is getting  work from shared queue, num dfs obtained = " << num_dfs << " regen level =" << gprv.embeddings_regeneration_level);

  for(int i = 0; i < gprv.embeddings_regeneration_level; i++) {
    if(gprv.dfs_task_queue.size() < (i + 1) ) {
      std::deque<types::DFS> tmp;
      //omp_set_lock(&lock);
      gprv.dfs_task_queue.push_back(tmp);
      //omp_unset_lock(&lock);
    }
    //omp_set_lock(&lock);
    types::DFS dfs = dfs_task_queue_shared[thread_id][i].back();
    gprv.dfs_task_queue[i].push_back(dfs);
    DEBUG(*(graph_miner::logger), "thread id " << thread_id << " inside thread process received data : " << dfs.to_string() << "  regen level " << gprv.embeddings_regeneration_level << ", queue size of level 0 = " << gprv.dfs_task_queue[0].size() );
    dfs_task_queue_shared[thread_id][i].pop_back();
    //omp_unset_lock(&lock);
  }

  if(gprv.dfs_task_queue.size() < ( gprv.embeddings_regeneration_level + 1) ) {
    std::deque<types::DFS> tmp;
    //omp_set_lock(&lock);
    gprv.dfs_task_queue.push_back(tmp);
    //omp_unset_lock(&lock);
  }

  for(int j = 0; j < num_dfs; j++) {
    //omp_set_lock(&lock);
    types::DFS dfs = dfs_task_queue_shared[thread_id][gprv.embeddings_regeneration_level].back();
    DEBUG(*(graph_miner::logger), "thread id " << thread_id << " inside thread process received data : " << dfs.to_string() << "  regen level " << gprv.embeddings_regeneration_level);
    gprv.dfs_task_queue[gprv.embeddings_regeneration_level].push_front(dfs);
    dfs_task_queue_shared[thread_id][gprv.embeddings_regeneration_level].pop_back();
    //omp_unset_lock(&lock);
  }

  //dfs_task_queue_shared[thread_id][gprv.embeddings_regeneration_level].clear();
  //thread_is_working[requesting_thread] = true;
  TRACE(*(graph_miner::logger), "thread " << thread_id << " embeddings regeneration level " << gprv.embeddings_regeneration_level);
}

void graph_miner_mpi_omp_hybrid::thread_start_working(){

  int thread_id = omp_get_thread_num();
  omp_set_lock(&lock);
  thread_is_working[thread_id] = true;
  DEBUG(*(graph_miner::logger), "thread " << thread_id << " is going to working status."  );
  omp_unset_lock(&lock);
}

void graph_miner_mpi_omp_hybrid::thread_split_global_work(int requester_rank_id, Thread_private_data &gprv){

  //half of the work to the global shared queue
  int thread_id = gprv.thread_id;       //omp_get_thread_num();
  DEBUG(*(graph_miner::logger), "inside thread split global work : " << gprv.DFS_CODE.to_string() << "  split level " << gprv.task_split_level << ", queue size of level 0 = " << gprv.dfs_task_queue[0].size() );

  if(!thread_working(gprv))
    return;

  gprv.task_split_level = 0;
  while(gprv.task_split_level < gprv.current_dfs_level && gprv.dfs_task_queue[gprv.task_split_level].size() < task_split_threshold ) {
    gprv.task_split_level++;
    DEBUG(*(graph_miner::logger), "thread" << thread_id << "task split level changed = " << gprv.task_split_level << " queue size = " << gprv.dfs_task_queue[gprv.task_split_level].size());
  }

  if( gprv.dfs_task_queue[gprv.task_split_level].size() < task_split_threshold )
    return;

  //DEBUG(*(graph_miner::logger), "thread" << thread_id << " final task split level = " <<gprv.task_split_level <<" donation size = "<< (gprv.dfs_task_queue[gprv.task_split_level].size()/2));

  for(int i = 0; i < gprv.task_split_level; i++) {
    if(global_shared_queue[requester_rank_id][thread_id].size() < (i + 1) ) {
      std::deque<types::DFS> tmp;
      global_shared_queue[requester_rank_id][thread_id].push_back(tmp);
    }
    global_shared_queue[requester_rank_id][thread_id][i].push_back(gprv.DFS_CODE[i]);
  }

  if(global_shared_queue[requester_rank_id][thread_id].size() < ( gprv.task_split_level + 1) ) {
    std::deque<types::DFS> tmp;
    global_shared_queue[requester_rank_id][thread_id].push_back(tmp);
  }

  int num_dfs = gprv.dfs_task_queue[gprv.task_split_level].size() / 2;

  DEBUG(*(graph_miner::logger), "thread " << thread_id << " is splitting work: split level = " << gprv.task_split_level << " queue size = " << gprv.dfs_task_queue[gprv.task_split_level].size() << " num dfs to donate = " << num_dfs << " requesting rank " << requester_rank_id);
  for(int j = 0; j < num_dfs; j++) {
    types::DFS dfs = gprv.dfs_task_queue[gprv.task_split_level].back();
    global_shared_queue[requester_rank_id][thread_id][gprv.task_split_level].push_front(dfs);
    DEBUG(*(graph_miner::logger), "thread id " << thread_id << " inside thread split work : " << dfs.to_string() << "  split level " << gprv.task_split_level << ", queue size of level 0 = " << gprv.dfs_task_queue[0].size() );
    gprv.dfs_task_queue[gprv.task_split_level].pop_back();
  }

}

void graph_miner_mpi_omp_hybrid::complete_global_work_split_request(int requester_rank_id, Thread_private_data &gprv){
  DEBUG(*(graph_miner::logger), "thread " << gprv.thread_id << " inside complete_global_work_split_request");
  complete_global_split_work(requester_rank_id, gprv);
}

} // namespace GRAPH_MINER

