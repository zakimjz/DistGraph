/*
   Copyright (C) 2004 Taku Kudo, All rights reserved.
   Copyright (C) 2014 Robert Kessl, All rights reserved.
   Copyright (C) 2015 Nilothpal Talukder, All rights reserved.
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

#include <iostream>
#include <dfs_code.hpp>
#include <graph_types.hpp>

using types::DFSCode;
using types::RMPath;
using types::Graph;
using types::Projected;
using types::History;

using types::Projected_map3;
using types::Projected_map2;
using types::Projected_map1;
using types::Projected_iterator3;
using types::Projected_iterator2;
using types::Projected_iterator1;


using std::cout;
using std::endl;

namespace types {

DFSCode Graph::get_min_dfs_code() const
{
  DFSCode min_dfs_code;

  Projected_map3 root;
  types::EdgeList edges;

  for(unsigned int from = 0; from < this->size(); ++from) {
    if(get_forward_root(*this, (*this)[from], edges)) {
      for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it) {
        root[(*this)[from].label][(*it)->elabel][(*this)[(*it)->to].label].push(0, *it, 0);
      } // for it
    } // if get_forward_root
  } // for from

  Projected_iterator3 fromlabel = root.begin();
  Projected_iterator2 elabel = fromlabel->second.begin();
  Projected_iterator1 tolabel = elabel->second.begin();

  min_dfs_code.push(0, 1, fromlabel->first, elabel->first, tolabel->first);

  get_min_dfs_code_internal(tolabel->second, min_dfs_code);
  return min_dfs_code;
} // dfs_code_is_min


void Graph::get_min_dfs_code_internal(Projected &projected, DFSCode &min_dfs_code) const
{
  const RMPath &rmpath = min_dfs_code.buildRMPath();
  int minlabel         = min_dfs_code[0].fromlabel;
  int maxtoc           = min_dfs_code[rmpath[0]].to;


  // SUBBLOCK 1
  {
    Projected_map1 root;
    bool flg = false;
    int newto = 0;

    // iterate over all backward edges that are connected to the right-most path
    // first found => proceed further in the dfs search for the minimum code.
    for(int i = rmpath.size() - 1; !flg  && i >= 1; --i) {
      for(unsigned int n = 0; n < projected.size(); ++n) {
        PDFS *cur = &projected[n];
        History history(*this, cur);
        Edge *e = get_backward(*this, history[rmpath[i]], history[rmpath[0]], history);
        if(e) {
          root[e->elabel].push(0, e, cur);
          newto = min_dfs_code[rmpath[i]].from;
          flg = true;
        } // if e
      } // for n
    } // for i

    // if we have found at least one backward edge, get the smallest and append its DFS to the min_dfs_code
    if(flg) {
      Projected_iterator1 elabel = root.begin();
      DFS new_min_dfs_code_elem(maxtoc, newto, -1, elabel->first, -1);
      min_dfs_code.push_back(new_min_dfs_code_elem);

      get_min_dfs_code_internal(elabel->second, min_dfs_code);
      return;
    } // if(flg)
  } // SUBBLOCK 1

  // SUBBLOCK 2
  {
    bool flg = false;
    int newfrom = 0;
    Projected_map2 root;
    types::EdgeList edges;

    // collect all forward edges
    // edges that extends the right-most node
    for(unsigned int n = 0; n < projected.size(); ++n) {
      PDFS *cur = &projected[n];
      History history(*this, cur);
      if(get_forward_pure(*this, history[rmpath[0]], minlabel, history, edges)) {
        flg = true;
        newfrom = maxtoc;
        for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
          root[(*it)->elabel][(*this)[(*it)->to].label].push(0, *it, cur);
      } // if get_forward_pure
    } // for n

    // edges that extends a node on the right-most path except the right-most node
    for(int i = 0; !flg && i < (int)rmpath.size(); ++i) {
      for(unsigned int n = 0; n < projected.size(); ++n) {
        PDFS *cur = &projected[n];
        History history(*this, cur);
        if(get_forward_rmpath(*this, history[rmpath[i]], minlabel, history, edges)) {
          flg = true;
          newfrom = min_dfs_code[rmpath[i]].from;
          for(types::EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
            root[(*it)->elabel][(*this)[(*it)->to].label].push(0, *it, cur);
        } // if get_forward_rmpath
      } // for n
    } // for i

    if(flg) {
      Projected_iterator2 elabel  = root.begin();
      Projected_iterator1 tolabel = elabel->second.begin();
      DFS new_min_dfs_code_elem(newfrom, maxtoc + 1, -1, elabel->first, tolabel->first);

      min_dfs_code.push_back(new_min_dfs_code_elem);
      get_min_dfs_code_internal(tolabel->second, min_dfs_code);

      return;
    } // if(flg)
  } // SUBBLOCK 2

  return;
} // Graph::get_min_dfs_code_internal


} // namespace types

