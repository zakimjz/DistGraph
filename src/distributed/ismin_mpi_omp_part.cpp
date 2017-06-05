/*
    $Id: ismin.cpp,v 1.5 2004/05/21 05:50:13 taku-ku Exp $;

   Copyright (C) 2004 Taku Kudo, All rights reserved.
   Copyright (C) 2015-2016 Nilothpal Talukder, All rights reserved.
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

#include <graph_miner_mpi_omp_part.hpp>
#include <types.hpp>

namespace GRAPH_MINER {

using types::Projected;
using types::Projected_map3;
using types::Projected_map2;
using types::Projected_map1;
using types::Projected_iterator3;
using types::Projected_iterator2;
using types::Projected_iterator1;
using types::EdgeList;
using types::History;
using types::PDFS;


bool graph_miner_mpi_omp_part::is_min(int thread_id)
{

  if(DFS_CODE_V[thread_id].size() == 1) {
    return (true);
  }

  DFS_CODE_V[thread_id].toGraph(GRAPH_IS_MIN_V[thread_id]);
  DFS_CODE_IS_MIN_V[thread_id].clear();

  Projected_map3 root;
  types::EdgeList edges;

  for(unsigned int from = 0; from < GRAPH_IS_MIN_V[thread_id].size(); ++from) {
    if(get_forward_root(GRAPH_IS_MIN_V[thread_id], GRAPH_IS_MIN_V[thread_id][from], edges)) {
      for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it) {
        root[GRAPH_IS_MIN_V[thread_id][from].label][(*it)->elabel][GRAPH_IS_MIN_V[thread_id][(*it)->to].label].push(0, *it, 0);
      } // for it
    } // if get_forward_root
  } // for from

  Projected_iterator3 fromlabel = root.begin();
  Projected_iterator2 elabel = fromlabel->second.begin();
  Projected_iterator1 tolabel = elabel->second.begin();

  DFS_CODE_IS_MIN_V[thread_id].push(0, 1, fromlabel->first, elabel->first, tolabel->first);

  return (project_is_min(thread_id,tolabel->second));
}

bool graph_miner_mpi_omp_part::project_is_min(int thread_id, Projected &projected)
{
  const RMPath& rmpath = DFS_CODE_IS_MIN_V[thread_id].buildRMPath();
  int minlabel         = DFS_CODE_IS_MIN_V[thread_id][0].fromlabel;
  int maxtoc           = DFS_CODE_IS_MIN_V[thread_id][rmpath[0]].to;


  // SUBBLOCK 1
  {
    Projected_map1 root;
    bool flg = false;
    int newto = 0;

    for(int i = rmpath.size() - 1; !flg  && i >= 1; --i) {
      for(unsigned int n = 0; n < projected.size(); ++n) {
        PDFS *cur = &projected[n];
        History history(GRAPH_IS_MIN_V[thread_id], cur);
        Edge *e = get_backward(GRAPH_IS_MIN_V[thread_id], history[rmpath[i]], history[rmpath[0]], history);
        if(e) {
          root[e->elabel].push(0, e, cur);
          newto = DFS_CODE_IS_MIN_V[thread_id][rmpath[i]].from;
          flg = true;
        } // if e
      } // for n
    } // for i

    if(flg) {
      Projected_iterator1 elabel = root.begin();
      DFS_CODE_IS_MIN_V[thread_id].push(maxtoc, newto, -1, elabel->first, -1);
      if(DFS_CODE_V[thread_id][DFS_CODE_IS_MIN_V[thread_id].size() - 1] != DFS_CODE_IS_MIN_V[thread_id][DFS_CODE_IS_MIN_V[thread_id].size() - 1]) return false;
      return project_is_min(thread_id, elabel->second);
    }
  } // SUBBLOCK 1

  // SUBBLOCK 2
  {
    bool flg = false;
    int newfrom = 0;
    Projected_map2 root;
    types::EdgeList edges;

    for(unsigned int n = 0; n < projected.size(); ++n) {
      PDFS *cur = &projected[n];
      History history(GRAPH_IS_MIN_V[thread_id], cur);
      if(get_forward_pure(GRAPH_IS_MIN_V[thread_id], history[rmpath[0]], minlabel, history, edges)) {
        flg = true;
        newfrom = maxtoc;
        for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
          root[(*it)->elabel][GRAPH_IS_MIN_V[thread_id][(*it)->to].label].push(0, *it, cur);
      } // if get_forward_pure
    } // for n

    for(int i = 0; !flg && i < (int)rmpath.size(); ++i) {
      for(unsigned int n = 0; n < projected.size(); ++n) {
        PDFS *cur = &projected[n];
        History history(GRAPH_IS_MIN_V[thread_id], cur);
        if(get_forward_rmpath(GRAPH_IS_MIN_V[thread_id], history[rmpath[i]], minlabel, history, edges)) {
          flg = true;
          newfrom = DFS_CODE_IS_MIN_V[thread_id][rmpath[i]].from;
          for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
            root[(*it)->elabel][GRAPH_IS_MIN_V[thread_id][(*it)->to].label].push(0, *it, cur);
        } // if get_forward_rmpath
      } // for n
    } // for i

    if(flg) {
      Projected_iterator2 elabel  = root.begin();
      Projected_iterator1 tolabel = elabel->second.begin();
      DFS_CODE_IS_MIN_V[thread_id].push(newfrom, maxtoc + 1, -1, elabel->first, tolabel->first);
      if(DFS_CODE_V[thread_id][DFS_CODE_IS_MIN_V[thread_id].size() - 1] != DFS_CODE_IS_MIN_V[thread_id][DFS_CODE_IS_MIN_V[thread_id].size() - 1]) return false;
      return project_is_min(thread_id, tolabel->second);
    } // if(flg)
  } // SUBBLOCK 2

  return true;
} // graph_miner::project_is_min


bool graph_miner_mpi_omp_part::is_min(Thread_private_data &gprv)
{

  if(gprv.DFS_CODE.size() == 1) {
    return (true);
  }

  gprv.DFS_CODE.toGraph(gprv.GRAPH_IS_MIN);
  gprv.DFS_CODE_IS_MIN.clear();

  Projected_map3 root;
  types::EdgeList edges;

  for(unsigned int from = 0; from < gprv.GRAPH_IS_MIN.size(); ++from) {
    if(get_forward_root(gprv.GRAPH_IS_MIN, gprv.GRAPH_IS_MIN[from], edges)) {
      for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it) {
        root[gprv.GRAPH_IS_MIN[from].label][(*it)->elabel][gprv.GRAPH_IS_MIN[(*it)->to].label].push(0, *it, 0);
      } // for it
    } // if get_forward_root
  } // for from

  Projected_iterator3 fromlabel = root.begin();
  Projected_iterator2 elabel = fromlabel->second.begin();
  Projected_iterator1 tolabel = elabel->second.begin();

  gprv.DFS_CODE_IS_MIN.push(0, 1, fromlabel->first, elabel->first, tolabel->first);

  return (project_is_min(gprv,tolabel->second));
}

bool graph_miner_mpi_omp_part::project_is_min(Thread_private_data &gprv, Projected &projected)
{
  const RMPath& rmpath = gprv.DFS_CODE_IS_MIN.buildRMPath();
  int minlabel         = gprv.DFS_CODE_IS_MIN[0].fromlabel;
  int maxtoc           = gprv.DFS_CODE_IS_MIN[rmpath[0]].to;


  // SUBBLOCK 1
  {
    Projected_map1 root;
    bool flg = false;
    int newto = 0;

    for(int i = rmpath.size() - 1; !flg  && i >= 1; --i) {
      for(unsigned int n = 0; n < projected.size(); ++n) {
        PDFS *cur = &projected[n];
        History history(gprv.GRAPH_IS_MIN, cur);
        Edge *e = get_backward(gprv.GRAPH_IS_MIN, history[rmpath[i]], history[rmpath[0]], history);
        if(e) {
          root[e->elabel].push(0, e, cur);
          newto = gprv.DFS_CODE_IS_MIN[rmpath[i]].from;
          flg = true;
        } // if e
      } // for n
    } // for i

    if(flg) {
      Projected_iterator1 elabel = root.begin();
      gprv.DFS_CODE_IS_MIN.push(maxtoc, newto, -1, elabel->first, -1);
      if(gprv.DFS_CODE[gprv.DFS_CODE_IS_MIN.size() - 1] != gprv.DFS_CODE_IS_MIN[gprv.DFS_CODE_IS_MIN.size() - 1]) return false;
      return project_is_min(gprv, elabel->second);
    }
  } // SUBBLOCK 1

  // SUBBLOCK 2
  {
    bool flg = false;
    int newfrom = 0;
    Projected_map2 root;
    types::EdgeList edges;

    for(unsigned int n = 0; n < projected.size(); ++n) {
      PDFS *cur = &projected[n];
      History history(gprv.GRAPH_IS_MIN, cur);
      if(get_forward_pure(gprv.GRAPH_IS_MIN, history[rmpath[0]], minlabel, history, edges)) {
        flg = true;
        newfrom = maxtoc;
        for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
          root[(*it)->elabel][gprv.GRAPH_IS_MIN[(*it)->to].label].push(0, *it, cur);
      } // if get_forward_pure
    } // for n

    for(int i = 0; !flg && i < (int)rmpath.size(); ++i) {
      for(unsigned int n = 0; n < projected.size(); ++n) {
        PDFS *cur = &projected[n];
        History history(gprv.GRAPH_IS_MIN, cur);
        if(get_forward_rmpath(gprv.GRAPH_IS_MIN, history[rmpath[i]], minlabel, history, edges)) {
          flg = true;
          newfrom = gprv.DFS_CODE_IS_MIN[rmpath[i]].from;
          for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
            root[(*it)->elabel][gprv.GRAPH_IS_MIN[(*it)->to].label].push(0, *it, cur);
        } // if get_forward_rmpath
      } // for n
    } // for i

    if(flg) {
      Projected_iterator2 elabel  = root.begin();
      Projected_iterator1 tolabel = elabel->second.begin();
      gprv.DFS_CODE_IS_MIN.push(newfrom, maxtoc + 1, -1, elabel->first, tolabel->first);
      if(gprv.DFS_CODE[gprv.DFS_CODE_IS_MIN.size() - 1] != gprv.DFS_CODE_IS_MIN[gprv.DFS_CODE_IS_MIN.size() - 1]) return false;
      return project_is_min(gprv, tolabel->second);
    } // if(flg)
  } // SUBBLOCK 2

  return true;
} // graph_miner::project_is_min


} // namespace GRAPH_MINER
