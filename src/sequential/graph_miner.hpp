/*
    $Id: gspan.h,v 1.6 2004/05/21 05:50:13 taku-ku Exp $;

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
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <types.hpp>
#include <graph_types.hpp>
#include <dfs_code.hpp>
#include <graph_output.hpp>
#include <logger.hpp>

namespace GRAPH_MINER {
using types::Edge;
using types::Vertex;
using types::Graph;
using types::DFS;
using types::DFSCode;
using types::RMPath;
using types::PDFS;
using types::History;
using types::Projected;

using types::Projected_map3;
using types::Projected_map2;
using types::Projected_map1;
using types::Projected_iterator3;
using types::Projected_iterator2;
using types::Projected_iterator1;

using types::Projected_riterator3;


class graph_miner {
protected:

  int minimal_support;
  graph_output * output;
  Logger *logger;

  types::Graph graph;
  DFSCode DFS_CODE;
  DFSCode DFS_CODE_IS_MIN;
  Graph GRAPH_IS_MIN;

  bool is_min();
  bool project_is_min(types::Projected &);

  std::map<unsigned int, unsigned int> support_counts(Projected &projected);
  virtual unsigned int support(Projected&);
  void project(Projected &);
  void report(Projected &, unsigned int);

  virtual void run_intern(void);

public:
  graph_miner();

  void set_graph(types::Graph &g);
  virtual void set_min_support(int minsup);
  virtual void set_graph_output(graph_output * gout);

  virtual void run();
};

} // namespace GRAPH_MINER


