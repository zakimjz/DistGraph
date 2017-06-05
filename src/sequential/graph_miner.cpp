/*
    $Id: gspan.cpp,v 1.8 2004/05/21 09:27:17 taku-ku Exp $;

   Copyright (C) 2004 Taku Kudo, All rights reserved.
   Copyright (C) 2014 Robert Kessl, All rights reserved.
   Copyright (C) 2014-2016 Nilothpal Talukder, All rights reserved.
     This is free software with ABSOLUTELY NO WARRANTY.

   This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.

   This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
     02111-1307, USA
 */
#include <graph_miner.hpp>
#include <iterator>

#include <stdlib.h>
#include <unistd.h>
#include <cstdio>
#include <iostream>

using namespace std;

namespace GRAPH_MINER {

graph_miner::graph_miner(void)
{
  logger = Logger::get_logger("GRAPH_MINER_SEQ");
  output = 0;
  minimal_support = -1;
}

void graph_miner::set_graph(types::Graph &g)
{
  this->graph = g;
}

void graph_miner::set_min_support(int minsup)
{
  minimal_support = minsup;
}

void graph_miner::set_graph_output(graph_output * gout)
{
  output = gout;
}


std::map<unsigned int, unsigned int>
graph_miner::support_counts(Projected &projected)
{
  std::map<unsigned int, unsigned int> counts;

  for(Projected::iterator cur = projected.begin();
      cur != projected.end(); ++cur) {
    counts[cur->id] += 1;
  }

  return (counts);
}

//support function for a single large graph, computes the minimum count of a node in the embeddings
unsigned int
graph_miner::support(Projected &projected)
{
  std::map<unsigned int, map<unsigned int, unsigned int> > node_id_counts;

  //Print DFS code
  //for(int i = 0; i < DFS_CODE.size(); i++)
  //  std::cout << DFS_CODE[i].to_string();
  //std::cout << endl;

  //iterated through the all the embeddings
  for(Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {
    PDFS *em = &(*cur);
    int dfsindex = DFS_CODE.size() - 1;
    while(em) {
      //std::cout<<DFS_CODE[dfsindex].to_string() <<endl;
      //print embedding
      //std::cout<<em->to_string()<<endl;

      if(DFS_CODE[dfsindex].to > DFS_CODE[dfsindex].from)    //forward edge
        node_id_counts[DFS_CODE[dfsindex].to][em->edge->to]++;
      if(!em->prev)
        node_id_counts[DFS_CODE[dfsindex].from][em->edge->from]++;
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
  //cout << "support: " << min << endl;
  return min;
}

void graph_miner::report(Projected &projected, unsigned int sup)
{
  output->output_graph(DFS_CODE, sup);
}

/* Recursive subgraph mining function (similar to subprocedure 1
 * Subgraph_Mining in [Yan2002]).
 */
void graph_miner::project(Projected &projected)
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

  DEBUG(*logger, "executing project for code: " << DFS_CODE.to_string() << "; support: " << sup);


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

  // Enumerate all possible one edge extensions of the current substructure.
  for(unsigned int n = 0; n < projected.size(); ++n) {

    unsigned int id = projected[n].id;
    PDFS *cur = &projected[n];
    History history(graph, cur);

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

  // Test all extended substructures.
  // backward
  for(Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
    for(Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
      DFS_CODE.push(maxtoc, to->first, -1, elabel->first, -1);
      project(elabel->second);
      DFS_CODE.pop();
    }
  }

  // forward
  for(Projected_riterator3 from = new_fwd_root.rbegin();
      from != new_fwd_root.rend(); ++from) {
    for(Projected_iterator2 elabel = from->second.begin();
        elabel != from->second.end(); ++elabel) {
      for(Projected_iterator1 tolabel = elabel->second.begin();
          tolabel != elabel->second.end(); ++tolabel) {
        DFS_CODE.push(from->first, maxtoc + 1, -1, elabel->first, tolabel->first);
        project(tolabel->second);
        DFS_CODE.pop();
      }
    }
  }

  return;
}

void graph_miner::run()
{
  run_intern();
}

void graph_miner::run_intern(void)
{

  types::EdgeList edges;
  Projected_map3 root;

  for(unsigned int from = 0; from < graph.size(); ++from) {
    if(get_forward_root(graph, graph[from], edges)) {   // get the edge list of the node g[from] in graph g
      for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
        root[graph[from].label][(*it)->elabel][graph[(*it)->to].label].push(0, *it, 0);    //projected (PDFS vector) entry: graph id (always 0 for single graph), edge pointer and null PDFS

    }   // if
  }   // for from
  //} // for id


  int total_single_edge_code = 0;

  for(Projected_iterator3 fromlabel = root.begin();
      fromlabel != root.end(); ++fromlabel) {
    for(Projected_iterator2 elabel = fromlabel->second.begin();
        elabel != fromlabel->second.end(); ++elabel) {
      for(Projected_iterator1 tolabel = elabel->second.begin();
          tolabel != elabel->second.end(); ++tolabel) {
        // Build the initial two-node graph.  It will be grownrecursively within project.
        DFS_CODE.push(0, 1, fromlabel->first, elabel->first, tolabel->first);
        project(tolabel->second); //Projected (PDFS vector): each entry contains graph id 0, edge pointer, null PDFS
        DFS_CODE.pop();
        total_single_edge_code++;
      } // for tolabel
    } // for elabel
  } // for fromlabel

  cout << " Total single edge DFS code " << total_single_edge_code << endl;
} // void graph_miner::run_intern(void)

} // namespace GRAPH_MINER


