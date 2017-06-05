/*
 * This source is part of the single graph mining algorithm.
 *
 * Copyright 2014 Robert Kessl
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

#ifndef __GRAPH_TYPES_HPP__
#define __GRAPH_TYPES_HPP__

#include <types.hpp>
#include <dfs_code.hpp>
#include <set>
#include <map>
#include <string>
#include <iostream>


namespace types {

typedef int edge_label_t;
typedef int vertex_label_t;

typedef std::set<edge_label_t> edge_label_set_t;
typedef std::set<vertex_label_t> vertex_label_set_t;

struct Edge {
  int from;
  int to;
  int elabel;
  unsigned int id;
  Edge() : from(0), to(0), elabel(0), id(0) {
  }
  std::string to_string() const;
  std::string to_string_global_vid(Graph &g) const;
};

class Vertex {
public:
  typedef std::vector<Edge>::iterator edge_iterator;
  typedef std::vector<Edge>::const_iterator const_edge_iterator;

  int label, global_vid, vertex_part_id, orig_part_id;
  bool is_boundary_vertex;
  std::vector<Edge> edge; //neighbor list

  void push(int from, int to, int elabel) {
    edge.resize(edge.size() + 1);
    edge[edge.size() - 1].from = from;
    edge[edge.size() - 1].to = to;
    edge[edge.size() - 1].elabel = elabel;
    return;
  }

  bool find(int from, int to, Edge &result) const {
    for(int i = 0; i < edge.size(); i++) {
      if(edge[i].from == from && edge[i].to == to) {
        result = edge[i];
        return true;
      }
    } // for i
    return false;
  } // find

  static size_t get_serialized_size(const Vertex &vrtx);
  static size_t get_serialized_size(char *buffer, size_t buffer_size);
  static size_t serialize(const Vertex &vrtx, char *buffer, size_t buffer_size);
  static size_t deserialize(Vertex &vrtx, char *buffer, size_t buffer_size);
};

class Graph : public std::vector<Vertex> {
private:
  unsigned int edge_size_;
public:
  typedef std::vector<Vertex>::iterator vertex_iterator;
  std::map<int,int> global_local_id_map;
  // for partitioned and external hop extension based approach this means the number of own vertices - 1
  // for all other approach this->size() - 1
  int max_local_vid;
  bool has_ext_neighbor;

  Graph(bool _directed);
  bool directed;


  //  int y; // class label
  unsigned int edge_size() const {
    return edge_size_;
  }
  unsigned int vertex_size() const {
    return (unsigned int)size();
  } // wrapper

  void buildEdge();
  std::istream &read_txt(std::istream &);  // read
  std::istream &read_adj(std::istream &);  // read
  std::ostream &write_txt(std::ostream &);  // write

  std::istream &read_partition_info(std::istream &, int part_id);  // read
  std::istream &read_partition_info_non_locals(std::istream &);  // read
  std::istream &read_adj_par(std::istream &);  // read

  std::istream &read_local_adjp(std::istream &is);
  std::ofstream &write_local_adjp(std::ofstream &os); // write

  void check(void);

  Graph();
  types::DFSCode get_min_dfs_code() const;
  int get_vertex_label(int vid) const {
    return at(vid).label;
  }

  int get_local_vid(int global_vid){
    if (global_local_id_map.count(global_vid) > 0)
      return global_local_id_map[global_vid];
    else
      return -1;
  }

  bool is_pseudo_local(int local_id){
    return (local_id > max_local_vid && (*this)[local_id].vertex_part_id != (*this)[local_id].orig_part_id);
  }

  bool is_external(int local_id){
    return (local_id > max_local_vid && (*this)[local_id].vertex_part_id == (*this)[local_id].orig_part_id);
  }

  bool is_subgraph(DFSCode &other_graph_dfs) const;
  void delete_edge(int from, int to);
  void delete_vertices(std::vector<int> local_ids);

  std::string to_string() const;

  static size_t get_serialized_size(const Graph &grph);
  static size_t get_serialized_size(char *buffer, size_t buffer_size);
  static size_t serialize(const Graph &grph, char *buffer, size_t buffer_size);
  static size_t deserialize(Graph &grph, char *buffer, size_t buffer_size);

  static size_t get_serialized_size(const graph_database_t &grph_db);
  static size_t get_serialized_size_db(char *buffer, size_t buffer_size);
  static size_t serialize(const graph_database_t &grph_db, char *buffer, size_t buffer_size);
  static size_t deserialize(graph_database_t &grph_db, char *buffer, size_t buffer_size);

  size_t get_serialized_size_for_partition(int vid);
  size_t serialize_neighbors_for_partition(int global_vid, int *buffer, int buffer_size);
  size_t get_serialized_size_for_partition(int vid, int requester_partition_id);
  size_t serialize_neighbors_for_partition(int requester_partition_id, int global_vid, int* &buffer, int buffer_size);
  size_t serialize_neighbors_for_partition(int requester_partition_id, int global_vid, int* &buffer, int buffer_size, std::set<int> &exclusions);
  size_t deserialize_neighbors_for_partition(int part_id, int*&buffer, int buffer_size, int& global_vid_recv);
  size_t deserialize_multiple_neighbors_for_partition(int part_id, int*&buffer, int buffer_size, int num_neighbors);

protected:
  void get_min_dfs_code_internal(Projected &projected, DFSCode &min_dfs_code) const;
  bool is_subgraph_internal(Projected &projected, DFSCode &other_graph_dfs, int depth) const;

};



struct PDFS {
  unsigned int id;      // ID of the original input graph
  Edge        *edge;
  PDFS        *prev;
  PDFS() : id(0), edge(0), prev(0) {
  };
  std::string to_string() const;
  //std::string to_string_projection(types::graph_database_t &gdb, types::graph_database_cuda &cgdb) const;
};


/**
 * Stores information of edges/nodes that were already visited in the
 * current DFS branch of the search.
 */

class History : public std::vector<Edge*> {
private:
  std::set<int> edge;
  std::set<int> vertex;

public:
  bool hasEdge(unsigned int id) {
    return (bool)edge.count(id);
  }
  bool hasVertex(unsigned int id) {
    return (bool)vertex.count(id);
  }
  void build(const Graph &, PDFS *);
  History() {
  }
  History(const Graph &g, PDFS *p) {
    build(g, p);
  }
  std::string to_string() const;
};


class Projected : public std::vector<PDFS> {
public:
  void push(int id, Edge *edge, PDFS *prev) {
    /*
        resize(size() + 1);
        //std::cout<< "capacity = " <<capacity()<< " size = "<<size() <<std::endl;
        //if(capacity() == (size() + 1) || capacity() == 0  )
        //  resize( (size() + 1) * 2 );
       PDFS &d = (*this)[size() - 1];
     */
    PDFS d;
    d.id = id;
    d.edge = edge;
    d.prev = prev;
    push_back(d);

  }
  std::string to_string() const;
};


struct Emb {
  Edge edge;
  Emb  *prev;
  std::string to_string() const;
  std::string to_string_global_vid(Graph &g);
};

class Embeddings : public std::vector<Emb> {
public:
  void push(Edge edge, Emb *prev) {
    Emb d;
    d.edge = edge;
    d.prev = prev;
    //if(prev)
    //std::cout<<" in push : prev "<< prev<< " "<<prev->to_string()<<std::endl;
    push_back(d);
  }
  std::string to_string() const;
  std::string print_global_vid(Graph &g);
};

struct Emb2 {
  Edge edge;
  int prev; //index into the previous vector of Emb2 or Embeddings2
  std::string to_string() const;
};

class Embeddings2 : public std::vector<Emb2> {
public:
  void push(Edge edge, int prev) {
    Emb2 d;
    d.edge = edge;
    d.prev = prev;
    //std::cout<<" in push : prev "<< prev<<std::endl;
    push_back(d);
  }
  std::string to_string() const;
};

class EmbVector : public std::vector<Edge> {
private:
  std::set<int> edge;
  std::set<int> vertex;

public:
  bool hasEdge(Edge e) {
    for(std::vector<Edge>::iterator it = this->begin(); it != this->end(); ++it) {
      if(it->from == e.from && it->to == e.to && it->elabel == e.elabel)
        return true;
      else if(it->from == e.to && it->to == e.from && it->elabel == e.elabel)
        return true;
    }
    return false;
    //return (bool)edge.count(id);
  }
  bool hasVertex(unsigned int id) {
    return (bool)vertex.count(id);
  }
  void build(const Graph &, Emb *);
  void build(const Graph &, const std::vector<Embeddings2>&, int, int);
  EmbVector() {
  }
  EmbVector(const Graph &g, Emb *p) {
    build(g, p);
  }
  EmbVector(const Graph &g, const std::vector<Embeddings2>& emb, int emb_col, int emb_index) {
    build(g, emb, emb_col, emb_index);
  }
  std::string to_string() const;
  std::string to_string_global_vid(Graph &g) const;
};


typedef std::map<int, std::map <int, std::map <int, Projected> > >           Projected_map3;
typedef std::map<int, std::map <int, Projected> >                            Projected_map2;
typedef std::map<int, Projected>                                             Projected_map1;
typedef std::map<int, std::map <int, std::map <int, Projected> > >::iterator Projected_iterator3;
typedef std::map<int, std::map <int, Projected> >::iterator Projected_iterator2;
typedef std::map<int, Projected>::iterator Projected_iterator1;
typedef std::map<int, std::map <int, std::map <int, Projected> > >::reverse_iterator Projected_riterator3;

typedef std::map<int, std::map <int, std::map <int, Embeddings> > >           Embeddings_map3;
typedef std::map<int, std::map <int, Embeddings> >                            Embeddings_map2;
typedef std::map<int, Embeddings>                                             Embeddings_map1;
typedef std::map<int, std::map <int, std::map <int, Embeddings> > >::iterator Embeddings_iterator3;
typedef std::map<int, std::map <int, Embeddings> >::iterator Embeddings_iterator2;
typedef std::map<int, Embeddings>::iterator Embeddings_iterator1;
typedef std::map<int, std::map <int, std::map <int, Embeddings> > >::reverse_iterator Embeddings_riterator3;

typedef std::map<int, std::map <int, std::map <int, Embeddings2> > >           Embeddings2_map3;
typedef std::map<int, std::map <int, Embeddings2> >                            Embeddings2_map2;
typedef std::map<int, Embeddings2>                                             Embeddings2_map1;
typedef std::map<int, std::map <int, std::map <int, Embeddings2> > >::iterator Embeddings2_iterator3;
typedef std::map<int, std::map <int, Embeddings2> >::iterator Embeddings2_iterator2;
typedef std::map<int, Embeddings2>::iterator Embeddings2_iterator1;
typedef std::map<int, std::map <int, std::map <int, Embeddings2> > >::reverse_iterator Embeddings2_riterator3;

typedef std::vector<int> graph_id_list_t;
typedef std::map<int, graph_id_list_t>   edge_gid_list1_t;
typedef std::map<int, edge_gid_list1_t>  edge_gid_list2_t;
typedef std::map<int, edge_gid_list2_t>  edge_gid_list3_t;


bool  get_forward_pure(const Graph &, Edge *,  int, History&, types::EdgeList &);
bool  get_forward_rmpath(const Graph &, Edge *,  int,  History&, types::EdgeList &);
bool  get_forward_root(const Graph &, const Vertex &, types::EdgeList &);
Edge *get_backward(const Graph &, Edge *,  Edge *, History&);

bool get_forward_pure(const Graph &, Edge,  int, EmbVector &, std::vector<Edge> &);
bool get_forward_rmpath(const Graph &, Edge,  int,  EmbVector &, std::vector<Edge> &);
bool  get_forward_root(const Graph &, const Vertex &, std::vector<Edge> &);
bool get_backward(const Graph &, Edge,  Edge, EmbVector &, Edge &);

bool  get_forward(const Graph &, const types::DFSCode &, History &, types::EdgeList &);
Edge *get_backward(const Graph &graph, const types::DFSCode &, History &);

bool  get_forward(const Graph &, const types::DFSCode &, EmbVector &, std::vector<Edge> &);
bool  get_backward(const Graph &graph, const types::DFSCode &, EmbVector &, Edge &);

} // namespace types

#endif

