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

#include <graph_miner_mpi_split_all.hpp>
#include <iterator>
#include <stdlib.h>
#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <math.h>

using namespace std;

namespace GRAPH_MINER {

graph_miner_mpi_split_all::graph_miner_mpi_split_all(void) : graph_miner()
{
  //std::cout<<"Inside graph_miner_mpi_dyn constructor"<<endl;
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  task_split_level = 0;   // splits are done from the top of the tree first
  task_split_threshold = 2;
  is_working = true;
  embeddings_regeneration_level = 0;
  current_dfs_level = 0;
  init_dynamic_load_balancing(rank, numtasks, MPI_COMM_WORLD);
}


void graph_miner_mpi_split_all::set_task_split_threshold(int thr){
  task_split_threshold = thr;
}

void graph_miner_mpi_split_all::set_load_balance_interval(int i){
  load_balance_interval = i;
}

/* Recursive subgraph mining function (similar to subprocedure 1
 * Subgraph_Mining in [Yan2002]).
 */
void graph_miner_mpi_split_all::project(Projected &projected, int dfs_level)
{
  // Check if the pattern is frequent enough.
  unsigned int sup = support(projected);

  if(sup < minimal_support) return;

  // The minimal DFS code check is more expensive than the support check,
  // hence it is done now, after checking the support.
  if(is_min() == false) {
    return;
  } else {

  }

  DEBUG(*logger, "DFS level = " << dfs_level);

  DEBUG(*(graph_miner::logger), "executing project for code: " << DFS_CODE.to_string() << "; support: " << sup);

  // Output the frequent substructure
  report(projected, sup);

  // In case we have a valid upper bound and our graph already exceeds it,
  // return.  Note: we do not check for equality as the DFS exploration may
  // still add edges within an existing subgraph, without increasing the
  // number of nodes.
  //
  //if(maxpat_max > maxpat_min && DFS_CODE.nodeCount() > maxpat_max) return;

  // We just outputted a frequent subgraph.  As it is frequent enough, so
  // might be its (n+1)-extension-graphs, hence we enumerate them all.
  const RMPath &rmpath = DFS_CODE.buildRMPath();
  int minlabel = DFS_CODE[0].fromlabel;
  int maxtoc = DFS_CODE[rmpath[0]].to;

  Projected_map3 new_fwd_root;
  Projected_map2 new_bck_root;
  types::EdgeList edges;

  current_dfs_level = dfs_level;

  // Enumerate all possible one edge extensions of the current substructure.
  for(unsigned int n = 0; n < projected.size(); ++n) {

    unsigned int id = projected[n].id;
    PDFS *cur = &projected[n];
    History history(graph, cur);

    // XXX: do we have to change something here for directed edges?

    // backward
    for(int i = (int)rmpath.size() - 1; i >= 1; --i) {
      Edge *e = get_backward(graph, history[rmpath[i]], history[rmpath[0]], history);
      if(e)
        new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push(id, e, cur);
    }

    // pure forward
    // FIXME: here we pass a too large e->to (== history[rmpath[0]]->to
    // into get_forward_pure, such that the assertion fails.
    //
    // The problem is:
    // history[rmpath[0]]->to > graph.size()
    if(get_forward_pure(graph, history[rmpath[0]], minlabel, history, edges)) {
      for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it) {
        new_fwd_root[maxtoc][(*it)->elabel][graph[(*it)->to].label].push(id, *it, cur);
      }
    }

    // backtracked forward
    for(int i = 0; i < (int)rmpath.size(); ++i) {
      if(get_forward_rmpath(graph, history[rmpath[i]], minlabel, history, edges)) {
        for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it) {
          new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][graph[(*it)->to].label].push(id, *it, cur);
        } // for it
      } // if
    } // for i
  } // for n


  std::deque<types::DFS> tmp;

  if(dfs_task_queue.size() <= dfs_level) {
    dfs_task_queue.push_back(tmp);
  }

  // Test all extended substructures.
  // backward
  for(Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
    for(Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {

      DFS dfs(maxtoc, to->first, -1, elabel->first, -1);
      dfs_task_queue[dfs_level].push_back(dfs);

      load_balance();
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
        dfs_task_queue[dfs_level].push_back(dfs);

        load_balance();
      }
    }
  }


  //current_dfs_level = dfs_level;
  //current_dfs_level = dfs_level + 1;

  while(dfs_task_queue[dfs_level].size() > 0) {

    DFS dfs = dfs_task_queue[dfs_level].front();
    dfs_task_queue[dfs_level].pop_front();
    DEBUG(*logger, "popped dfs = " << dfs.to_string() );

    current_dfs_level = dfs_level;
    load_balance();

    DFS_CODE.push(dfs.from, dfs.to, dfs.fromlabel, dfs.elabel, dfs.tolabel);

    if(dfs.is_backward())
      project(new_bck_root[dfs.to][dfs.elabel], dfs_level + 1);      //Projected (PDFS vector): each entry contains graph id 0, edge pointer, null PDFS
    else
      project(new_fwd_root[dfs.from][dfs.elabel][dfs.tolabel], dfs_level + 1);      //Projected (PDFS vector): each entry contains graph id 0, edge pointer, null PDFS

    DFS_CODE.pop();
  }

  //current_dfs_level = dfs_level;

  return;
}

/* Regnerate the embeddings for the dfs codes obtained from other processors
 */

void graph_miner_mpi_split_all::regenerate_embeddings(Projected &projected, int dfs_level)
{
  // We don't need to check if the pattern is frequent or minimal

  DEBUG(*(graph_miner::logger), "DFS level inside regenerate embeddings = " << dfs_level << " queue size = " << dfs_task_queue[dfs_level].size());

  //not necessary though, as task split is not done while regenerating embeddings

  //current_dfs_level = dfs_level + 1;

  //iterate for all in the task_queue

  for(int i = 0; dfs_task_queue[dfs_level].size() > 0; i++) {


    types::DFS dfs = dfs_task_queue[dfs_level].front();
    dfs_task_queue[dfs_level].pop_front();

    current_dfs_level = dfs_level;
    load_balance();

    DFS_CODE.push(dfs.from, dfs.to, dfs.fromlabel, dfs.elabel, dfs.tolabel);

    DEBUG(*(graph_miner::logger), "*****regenerating embeddings for code: " << DFS_CODE.to_string() );

    //const RMPath &rmpath = DFS_CODE.buildRMPath();
    //int minlabel = DFS_CODE[0].fromlabel;
    //int maxtoc = DFS_CODE[rmpath[0]].to;

    Projected new_root;

    for(unsigned int n = 0; n < projected.size(); ++n) {

      unsigned int id = projected[n].id;
      PDFS *cur = &projected[n];
      History history(graph, cur);

      if(dfs.is_backward() ) {
        Edge *e = get_backward(graph, DFS_CODE, history);
        if(e)
          new_root.push(id, e, cur);
      }else{
        types::EdgeList edges;
        if(get_forward(graph, DFS_CODE, history, edges)) {
          for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it) {
            new_root.push(id, *it, cur);
          }
        }
      }
    }

    if( embeddings_regeneration_level > dfs_level ) {
      regenerate_embeddings(new_root, dfs_level + 1);
    }else{
      //regeneration of embeddings ended
      //now perform regular extensions with project function
      //reset embeddings_regeneration_level
      //embeddings_regeneration_level = 0;
      project(new_root, dfs_level + 1);
    }

    DFS_CODE.pop();

  }

  //New in "split task in all levels"
  // decrement the embeddings_regeration_level after done with the queue
  embeddings_regeneration_level--;

  return;
}


void graph_miner_mpi_split_all::run_intern(void)
{

  types::EdgeList edges;
  Projected_map3 root;
  int single_edge_dfscodes = 0;

  // id is always 0 for a single graph

  //for(unsigned int id = 0; id < graph_database.size(); ++id) {
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


  int dfscodes_per_rank =  (int) ceil((single_edge_dfscodes * 1.0) / numtasks);
  int start_index = rank * dfscodes_per_rank;
  int end_index = start_index + dfscodes_per_rank - 1;
  if (end_index > single_edge_dfscodes - 1)
    end_index = single_edge_dfscodes - 1;

  DEBUG(*(graph_miner::logger), "start index = " << start_index << " , end index = " << end_index << endl);

  std::deque<types::DFS> tmp;
  dfs_task_queue.push_back(tmp);

  int index = 0;
  for(Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end(); ++fromlabel) {
    for(Projected_iterator2 elabel = fromlabel->second.begin(); elabel != fromlabel->second.end(); ++elabel) {
      for(Projected_iterator1 tolabel = elabel->second.begin();
          tolabel != elabel->second.end(); ++tolabel) {

        if( index >= start_index && index <= end_index ) {
          // Build the initial two-node graph.  It will be grownrecursively within project.

          DFS dfs(0, 1, fromlabel->first, elabel->first, tolabel->first);
          dfs_task_queue[0].push_back(dfs);
          //std::cout << dfs.to_string() << endl;
        }
        index++;

      } // for tolabel
    } // for elabel
  } // for fromlabel

  //std::cout<<"size = " << dfs_task_queue[0].size() << std::endl;

  //while(dfs_task_queue[0].size() > 0){
  while(computation_end == false) {

    if(dfs_task_queue[0].size() == 0) {
      is_working = false;
      embeddings_regeneration_level = 0;
      task_split_level = 0;
      load_balance();

    }else{
      //this is done in process_received_data, so not required here
      //is_working = true;

      DFS dfs = dfs_task_queue[0].front();
      dfs_task_queue[0].pop_front();
      DEBUG(*(graph_miner::logger), "popped dfs = " << dfs.to_string() );
      load_balance();

      DFS_CODE.push(0, 1, dfs.fromlabel, dfs.elabel, dfs.tolabel);
      current_dfs_level = 1;

      //INFO(*(graph_miner::logger), "embeddings regeneration level = " << embeddings_regeneration_level);
      if(embeddings_regeneration_level < 1)
        project(root[dfs.fromlabel][dfs.elabel][dfs.tolabel], 1);                    //Projected (PDFS vector): each entry contains graph id 0, edge pointer, null PDFS
      else
        regenerate_embeddings(root[dfs.fromlabel][dfs.elabel][dfs.tolabel], 1);

      current_dfs_level = 0;
      DFS_CODE.pop();
      if(dfs_task_queue[0].size() == 0) {
        DEBUG(*(graph_miner::logger),"processor " << rank << " is idle, has token = " << has_token);
        if(has_token == true)
          DEBUG(*(graph_miner::logger),"processor " << rank << " token color = " << get_token_color(token_color));
      }
    }

  }

} // void graph_miner_mpi_dyn::run_intern(void)

/* Implemented functions of dynamic_load_balancing class */

bool graph_miner_mpi_split_all::working(){
  return is_working;
}

bool graph_miner_mpi_split_all::can_split_work(){

  //This will prevent splitting work beyond the first level (single edge)
  //if(task_split_level == 0 && dfs_task_queue[task_split_level].size() == 0)
  //return false;

  if(!working())
    return false;

  //prevents giving work until regeneration is complete
  //if(embeddings_regeneration_level > 0)
  //	return false;

  task_split_level = current_dfs_level;
  //while(dfs_task_queue[task_split_level].size() < task_split_threshold && task_split_level < dfs_task_queue.size() ){
  while(dfs_task_queue[task_split_level].size() < task_split_threshold && task_split_level > 0 ) {
    task_split_level--;
    DEBUG(*(graph_miner::logger), "task split level changed = " << task_split_level << " queue size = " << dfs_task_queue[task_split_level].size());
  }


  DEBUG(*(graph_miner::logger), "task split level = " << task_split_level << " queue sizes = ");


  if(dfs_task_queue[task_split_level].size() >= task_split_threshold)
    return true;

  return false;
}


void graph_miner_mpi_split_all::split_work( int* &buffer, int &length){

  //send half of the queue for each level starting from root (level 0) to task_split_level

  //init send count for each level
  vector<int> num_dfs_for_levels;
  for(int i = 0; i <= task_split_level; i++)
    num_dfs_for_levels.push_back(0);

  //total count
  int num_dfs_to_send = 0;

  //consider the prefix code and split the task queue associated with these levels
  for(int i = 0; i < task_split_level; i++) {
    //consider the currently expanded i-th level prefix already in DFS_CODE
    num_dfs_for_levels[i] += 1;
    num_dfs_to_send += 1;
    //check if splittable and split
    if(dfs_task_queue[i].size() / 2 >= task_split_threshold) {
      num_dfs_for_levels[i] += dfs_task_queue[i].size() / 2;
      num_dfs_to_send += dfs_task_queue[i].size() / 2;
    }
  }

  //now split the queue of task_split_level (prefix not required)
  num_dfs_for_levels[task_split_level] = dfs_task_queue[task_split_level].size() / 2;
  num_dfs_to_send += dfs_task_queue[task_split_level].size() / 2;

  // allocate memory for storing task_split_level, counts in each level, and dfs codes in the queue
  buffer = new int[ 1 + (task_split_level + 1) + 5 * num_dfs_to_send];
  buffer[0] = task_split_level;
  std::stringstream ss1(" ");

  for(int i = 0; i <= task_split_level; i++) {
    buffer[1 + i] = num_dfs_for_levels[i];
    ss1 << buffer[1 + i] << " ";
  }

  DEBUG(*(graph_miner::logger), " sending dfs codes per level: " << ss1.str());
  //compute the initial offset for inserting dfs codes into buffer
  int offset = 1 + (task_split_level + 1);

  std::string str1 = "";
  int i = 0, j = 0;
  for(i = 0; i <= task_split_level; i++) {
    //add the trailing ones from the task queue first, ignore the prefix (unrolling will be done in reverse)
    for(j = 0; j < num_dfs_for_levels[i] - 1; j++) {

      types::DFS dfs = dfs_task_queue[i].back();
      DEBUG(*(graph_miner::logger), "DFS code to be sent " << dfs.to_string());
      dfs_task_queue[i].pop_back();
      dfs.serialize( (char*) (buffer + offset + j * 5), sizeof(int) * 5 );
      //str2 += str1 + dfs.to_string()+ ", ";

      //if(task_split_level > 0)
      //		str1 += DFS_CODE[i].to_string();
    }

    if( i < task_split_level) {
      //add the prefix code at the back
      DFS_CODE[i].serialize( (char*) ((buffer + offset) + j * 5), sizeof(int) * 5 );
    }else{
      //no prefix for the deepest level, so add one more from the task queue
      types::DFS dfs = dfs_task_queue[i].back();
      DEBUG(*(graph_miner::logger), "DFS code to be sent " << dfs.to_string());
      dfs_task_queue[i].pop_back();
      dfs.serialize( (char*) (buffer + offset + j * 5), sizeof(int) * 5 );
      //str2 += str1 + dfs.to_string()+ ", ";
    }

    //compute the new offset
    offset += num_dfs_for_levels[i] * 5;

  }

  DEBUG(*(graph_miner::logger), "Sent DFS code list" << str2);

  length = 1 + (task_split_level + 1) + 5 * num_dfs_to_send;

  //DEBUG(*(graph_miner::logger), "The size of dfs task queue after split = " << dfs_task_queue[0].size() );

}

void graph_miner_mpi_split_all::process_received_data(int* buffer, int size){

  int received_dfs_level = buffer[0];
  vector<int> num_dfs_for_levels;
  int num_dfs_received = 0;

  std::stringstream ss1("");
  for(int i = 0; i <= received_dfs_level; i++) {
    num_dfs_for_levels.push_back(buffer[1 + i]);
    num_dfs_received += buffer[1 + i];
    ss1 << buffer[1 + i] << " ";
  }

  DEBUG(*(graph_miner::logger), "Processing received data of size " << size);
  DEBUG(*(graph_miner::logger), "received dfs level = " << received_dfs_level << " num dfs codes = " << num_dfs_received );
  std::stringstream ss("Received data: ");
  DEBUG(*(graph_miner::logger), " received dfs codes per level: " << ss1.str());

  int offset = 1 + (received_dfs_level  + 1);
  int i = 0, j = 0;
  for(i = 0; i <= received_dfs_level; i++) {

    if(dfs_task_queue.size() < (i + 1) ) {
      std::deque<types::DFS> tmp;
      dfs_task_queue.push_back(tmp);
    }

    for(j = 0; j < num_dfs_for_levels[i]; j++) {
      types::DFS dfs;
      dfs.deserialize((char*) (buffer + offset + j * 5), 5 * sizeof(int));
      dfs_task_queue[i].push_front(dfs);
      ss << " level = " << i << " " << dfs_task_queue[i].front() << ", ";
      DEBUG( *(graph_miner::logger), "level = " << i << " queue size = " << dfs_task_queue[i].size());
    }
    //compute the new offset
    offset += num_dfs_for_levels[i] * 5;
    DEBUG(*(graph_miner::logger), " new offset : " << offset);
  }

  //DEBUG( *(graph_miner::logger), "level = " << received_dfs_level <<" queue size = "<<dfs_task_queue[received_dfs_level].size());
  //this check is not required
  //if(dfs_task_queue[received_dfs_level].size() > 0){
  is_working = true;
  DEBUG(*(graph_miner::logger), "processor " << rank << "is going to working status. has token = " << has_token);
  //}

  embeddings_regeneration_level = received_dfs_level;
  DEBUG(*(graph_miner::logger), ss.str());

  delete[] buffer;
}

} // namespace GRAPH_MINER

