/*
/*
 * This source is part of the single graph mining algorithm.
 * Copyright 2004 Taku Kudo
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
#include <graph_types.hpp>
#include <algorithm>
#include <cassert>
#include <string>
#include <sstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <logger.hpp>
#include <utils.hpp>
#include <cstring>
#include <mpi.h>

namespace types {

template <class T, class Iterator>
void tokenize(const char *str, Iterator iterator)
{
  std::istringstream is(std::string(str));
  std::copy(std::istream_iterator <T> (is), std::istream_iterator <T> (), iterator);
}

Graph::Graph() : edge_size_(0), directed(false)
{

}

Graph::Graph(bool _directed)
{
  directed = _directed;
}

void Graph::buildEdge()
{
  char buf[512];
  std::map <std::string, unsigned int> tmp;

  unsigned int id = 0;
  for(int from = 0; from < (int)size(); ++from) {
    for(Vertex::edge_iterator it = (*this)[from].edge.begin();
        it != (*this)[from].edge.end(); ++it) {
      if(directed || from <= it->to)
        std::sprintf(buf, "%d %d %d", from, it->to, it->elabel);
      else
        std::sprintf(buf, "%d %d %d", it->to, from, it->elabel);

      // Assign unique id's for the edges.
      if(tmp.find(buf) == tmp.end()) {
        it->id = id;
        tmp[buf] = id;
        ++id;
      } else {
        it->id = tmp[buf];
      }
    }
  }

  edge_size_ = id;
}


std::istream &Graph::read_txt(std::istream &is)
{
  char line[1024];
  std::vector<std::string> result;

  clear();

  while(true) {
    unsigned int pos = is.tellg();
    if(!is.getline(line, 1024)) {
      break;
    }
    result.clear();
    utils::split(line, result);

    if(result.empty()) {
      // do nothing
    } else if(result[0] == "t") {
      if(!empty()) {   // use as delimiter
        is.seekg(pos, std::ios_base::beg);
        break;
      } else {
        // y = atoi (result[3].c_str());
      }
    } else if(result[0] == "v" && result.size() >= 3) {
      unsigned int id    = atoi(result[1].c_str());
      this->resize(id + 1);
      (*this)[id].label = atoi(result[2].c_str());
      (*this)[id].global_vid = id;
    } else if(result[0] == "e" && result.size() >= 4) {
      int from   = atoi(result[1].c_str());
      int to     = atoi(result[2].c_str());
      int elabel = atoi(result[3].c_str());

      if((int)size() <= from || (int)size() <= to) {
        std::cerr << "Format Error:  define vertex lists before edges, from: " << from << "; to: " << to << "; vertex count: " << size() << std::endl;
        throw std::runtime_error("Format Error:  define vertex lists before edges");
      }

      (*this)[from].push(from, to, elabel);
      if(directed == false)
        (*this)[to].push(to, from, elabel);
    }
  }

  this->max_local_vid = this->size() - 1;
  buildEdge();

  return is;
}


std::ostream &Graph::write_txt(std::ostream &os)
{
  char buf[512];
  std::set <std::string> tmp;

  for(int from = 0; from < (int)size(); ++from) {
    os << "v " << from << " " << (*this)[from].label << std::endl;

    for(Vertex::edge_iterator it = (*this)[from].edge.begin();
        it != (*this)[from].edge.end(); ++it) {
      if(directed || from <= it->to) {
        std::sprintf(buf, "%d %d %d", from, it->to,   it->elabel);
      } else {
        std::sprintf(buf, "%d %d %d", it->to,   from, it->elabel);
      }
      tmp.insert(buf);
    } // for it
  } // for from

  for(std::set<std::string>::iterator it = tmp.begin(); it != tmp.end(); ++it) {
    os << "e " << *it << std::endl;
  } // for it

  return os;
}


std::ofstream &Graph::write_local_adjp(std::ofstream &os)
{
  os << (*this).size() << std::endl;
  std::cout << " total vertices = " << (*this).size() << std::endl;

  for(int from = 0; from < (int)size(); ++from) {
    os <<  (*this)[from].vertex_part_id << " " << (*this)[from].global_vid << " " << (*this)[from].label;

    for(Vertex::edge_iterator it = (*this)[from].edge.begin();
        it != (*this)[from].edge.end(); ++it) {
      os << " " << it->to  << " " << it->elabel;
    } // for it
    os << std::endl;
  } // for from

  return os;
}

/*
   This reads the adjacency list format also used by METIS graph partitioner
   Note: The vertex_ids and labels both start from 1, instead of 0
 */
std::istream &Graph::read_adj(std::istream &is)
{
  char line[1024];
  std::vector<std::string> result;

  clear();

  int num_vertices, num_edges, vertex_id;
  unsigned int pos = is.tellg();
  is.getline(line, 1024);
  result.clear();
  utils::split(line, result);
  if(result.empty()) {
    std::cerr << "Empty first line" << endl;
  } else {
    num_vertices = atoi(result[0].c_str());
    num_edges = atoi(result[1].c_str());
  }


  vertex_id = 0;
  while(true) {
    pos = is.tellg();
    if(!is.getline(line, 1024)) {
      break;
    }
    result.clear();
    utils::split(line, result);

    if(result.empty()) {
      // do nothing
    } else {
      this->resize(vertex_id + 1);
      (*this)[vertex_id].label = atoi(result[0].c_str()) - 1;
      (*this)[vertex_id].global_vid = vertex_id;
      for(int i = 1; i < result.size(); i++) {
        int to = atoi(result[i++].c_str()) - 1;
        int elabel = atoi(result[i].c_str()) - 1;
        (*this)[vertex_id].push(vertex_id, to, elabel);

      }

    }
    vertex_id++;
  }

  this->max_local_vid = this->size() - 1;
  buildEdge();

  return is;
}


/*
   This reads the adjacency list format also used by METIS graph partitioner
   Note: The vertex_ids and labels both start from 1, instead of 0
 */
std::istream &Graph::read_local_adjp(std::istream &is)
{
  char line[10240];
  std::vector<std::string> result;

  clear();

  int num_vertices, num_edges, vertex_id;
  unsigned int pos = is.tellg();
  is.getline(line, 10240);
  result.clear();
  utils::split(line, result);
  if(result.empty()) {
    std::cerr << "Empty first line" << endl;
  } else {
    num_vertices = atoi(result[0].c_str());
    //std::cout << " number of vertices " <<num_vertices<<endl;
  }

  this->resize(num_vertices); // allocate max possible
  vertex_id = 0;
  int part_id;
  this->max_local_vid = -1;
  while(true) {
    pos = is.tellg();
    if(!is.getline(line, 10240)) {
      break;
    }
    result.clear();
    utils::split(line, result);

    if(result.empty()) {
      // do nothing
    } else {
      types::Vertex v;
      this->push_back(v);
      //this->resize(vertex_id + 1);
      (*this)[vertex_id].vertex_part_id = atoi(result[0].c_str());
      (*this)[vertex_id].orig_part_id = (*this)[vertex_id].vertex_part_id;
      (*this)[vertex_id].global_vid = atoi(result[1].c_str());
      (*this)[vertex_id].label = atoi(result[2].c_str());
      (*this).global_local_id_map[(*this)[vertex_id].global_vid] = vertex_id;
      for(int i = 3; i < result.size(); i++) {
        int to = atoi(result[i++].c_str());
        int elabel = atoi(result[i].c_str());
        (*this)[vertex_id].push(vertex_id, to, elabel);
      }

    }

    //TODO: provide it inside the input file
    if(vertex_id == 0)
      part_id = (*this)[vertex_id].vertex_part_id;
    if((*this)[vertex_id].vertex_part_id == part_id)
      max_local_vid++;
    else if (this->has_ext_neighbor == false)
      this->has_ext_neighbor = true;

    vertex_id++;
  }
  //std::cout << " number of vertices " <<num_vertices<< " max local id " << max_local_vid << endl;
  this->resize(vertex_id);
  buildEdge();

  return is;
}


std::istream &Graph::read_adj_par(std::istream &is)
{
  char line[10240];
  std::vector<std::string> result;
  std::vector<int> non_locals_list;

  int num_vertices, num_edges, vertex_id;
  unsigned int pos = is.tellg();
  is.getline(line, 10240);
  result.clear();
  utils::split(line, result);
  if(result.empty()) {
    std::cerr << "Empty first line" << endl;
  } else {
    num_vertices = atoi(result[0].c_str());
    num_edges = atoi(result[1].c_str());
  }

  cout << "total vertices = " << num_vertices << " total edges = " << num_edges << endl;

  int max_local_id = this->size() - 1;
  this->max_local_vid = max_local_id;
  int global_vid = 0, local_vid = 0;
  int non_locals = 0;
  int int_edges = 0, ext_edges = 0;

  while(local_vid <= max_local_id) {

    pos = is.tellg();
    if(!is.getline(line, 10240)) {
      break;
    }
    result.clear();
    utils::split(line, result);

    if (global_vid == (*this)[local_vid].global_vid) {
      if(result.empty()) {
        // do nothing
      } else {

        (*this)[local_vid].label = atoi(result[0].c_str()) - 1;

        for(int i = 1; i < result.size(); i++) {

          int global_to = atoi(result[i++].c_str()) - 1;
          int elabel = atoi(result[i].c_str()) - 1;

          if( (*this).global_local_id_map.count(global_to) == 0) {

            non_locals_list.push_back(global_to);
            this->resize(max_local_id +  non_locals + 2);
            // we also need to get the tolabel
            (*this)[max_local_id + non_locals + 1].label = 0; // wrong label, fix it in later pass
            (*this)[max_local_id + non_locals + 1].global_vid = global_to;
            (*this).global_local_id_map[global_to] = max_local_id + non_locals + 1;

            non_locals++;
          }

          int to = (*this).global_local_id_map[global_to];
          (*this)[local_vid].push(local_vid, to, elabel);
          if( to > max_local_id) {
            (*this)[to].push(to, local_vid, elabel);
            ext_edges++;
          }else{
            int_edges++;
          }
        }
      }
      local_vid++;
    }
    global_vid++;
  }

  //cout << "read adj par: max local id = " << local_vid << ",  max global id = " << global_vid << ", non local vertices = " << non_locals << ", internal edges = "<< (int_edges/2) << ", external edges = " << ext_edges << endl;
  cout << "read adj par: internal vertices = " << local_vid << ", external vertices = " << non_locals << ", internal edges = " << (int_edges / 2) << ", external edges = " << ext_edges << endl;
  //read the graph file again to get the non-local vertex labels
  std::sort(non_locals_list.begin(), non_locals_list.end());

  is.clear();
  is.seekg(0,std::ios_base::beg);

  pos = is.tellg();
  is.getline(line, 1024);
  //skip first line

  global_vid = 0;
  int index = 0;

  while(index < non_locals_list.size()) {

    pos = is.tellg();
    if(!is.getline(line, 1024)) {
      break;
    }

    result.clear();
    utils::split(line, result);

    if (global_vid == non_locals_list[index]) {
      if(result.empty()) {
        // do nothing
      } else {
        local_vid = (*this).global_local_id_map[global_vid];
        (*this)[local_vid].label =  atoi(result[0].c_str()) - 1;
      }
      index++;
    }
    global_vid++;
  }

  buildEdge();

  return is;

}

std::istream &Graph::read_partition_info(std::istream &is, int part_id)
{
  char line[1024];
  std::vector<std::string> result;

  clear();

  int global_vid = 0, local_vid = 0;
  while(true) {
    unsigned int pos = is.tellg();
    if(!is.getline(line, 1024)) {
      break;
    }
    result.clear();
    utils::split(line, result);

    if(result.empty()) {
      // do nothing
    } else if (part_id ==  atoi(result[0].c_str())) {
      this->resize(local_vid + 1);
      (*this)[local_vid].global_vid = global_vid;
      //  cout << "read partition info: global id " << global_vid << " partition " << part_id <<endl;
      (*this)[local_vid].vertex_part_id = part_id;
      (*this)[local_vid].orig_part_id = part_id;
      (*this)[local_vid].is_boundary_vertex = false;
      this->global_local_id_map[global_vid] = local_vid;
      local_vid++;
    } else if (this->has_ext_neighbor == false) {
      this->has_ext_neighbor = true;
    }
    global_vid++;
  }

  cout << this->size() << " " << local_vid << " " << global_vid << endl;
  //buildEdge();

  this->max_local_vid = this->size() - 1;

  return is;
}

std::istream &Graph::read_partition_info_non_locals(std::istream &is)
{
  char line[1024];
  std::vector<std::string> result;

  int global_vid = 0, local_vid = 0;
  while(true) {
    unsigned int pos = is.tellg();
    if(!is.getline(line, 1024)) {
      break;
    }
    result.clear();
    utils::split(line, result);

    if(result.empty()) {
      // do nothing
    } else if (global_local_id_map.count(global_vid) > 0 ) {
      int part_id =  atoi(result[0].c_str());
      local_vid = global_local_id_map[global_vid];
      (*this)[local_vid].vertex_part_id = part_id;
      (*this)[local_vid].orig_part_id = part_id;
      //std::cout << "global id " << global_vid << " local id " << local_vid << " partition " <<part_id<<endl;
    }
    global_vid++;
  }

  return is;
}

void Graph::check(void)
{
  // Check all indices
  for(int from = 0; from < (int)size(); ++from) {
    //mexPrintf ("check vertex %d, label %d\n", from, (*this)[from].label);

    for(Vertex::edge_iterator it = (*this)[from].edge.begin();
        it != (*this)[from].edge.end(); ++it) {
      //mexPrintf ("   check edge from %d to %d, label %d\n", it->from, it->to, it->elabel);
      assert(it->from >= 0 && it->from < size());
      assert(it->to >= 0 && it->to < size());
    }
  }
}

std::string Edge::to_string() const {
  std::stringstream ss;
  ss << "e(" << from << "," << to << "," << elabel << ")";
  return ss.str();
}

std::string Edge::to_string_global_vid(Graph &g) const {
  std::stringstream ss;
  ss << "e(" << g[from].global_vid << "," << g[to].global_vid << "," << elabel << ")";
  return ss.str();
}

std::string PDFS::to_string() const {
  std::stringstream ss;
  ss << "[" << id << "," << edge->to_string() << "]";
  return ss.str();
}

std::string Projected::to_string() const
{
  std::stringstream ss;

  for(int i = 0; i < size(); i++) {
    ss << (*this)[i].to_string() << "; ";
  } // for i

  return ss.str();
} // Projected::to_string

std::string Emb::to_string() const {
  std::stringstream ss;
  ss << "[" << edge.to_string() << " " << prev << "]";
  return ss.str();
}

std::string Emb::to_string_global_vid(Graph &g) {
  Emb *emb = this;
  EmbVector ev = EmbVector(g, emb);
  return ev.to_string_global_vid(g);
}

std::string Emb2::to_string() const {
  std::stringstream ss;
  ss << "[" << edge.to_string() << " " << prev << "]";
  return ss.str();
}

std::string Embeddings::to_string() const
{
  std::stringstream ss;

  for(int i = 0; i < size(); i++) {
    ss << (*this)[i].to_string() << "; ";
  } // for i

  return ss.str();
} // Embeddings::to_string

std::string Embeddings::print_global_vid(Graph &g) {
  for(std::vector<Emb>::iterator it = this->begin(); it != this->end(); ++it) {
    std::cout << (*it).to_string_global_vid(g) << std::endl;
  }
}

std::string Embeddings2::to_string() const
{
  std::stringstream ss;

  for(int i = 0; i < size(); i++) {
    ss << (*this)[i].to_string() << "; ";
  } // for i

  return ss.str();
} // Embeddings::to_string


void History::build(const Graph &graph, PDFS *e)
{
  // first build history

  if(e) {
    push_back(e->edge);
    edge.insert(e->edge->id);
    vertex.insert(e->edge->from);
    vertex.insert(e->edge->to);

    for(PDFS *p = e->prev; p; p = p->prev) {
      push_back(p->edge);       // this line eats 8% of overall instructions(!)
      edge.insert(p->edge->id);
      vertex.insert(p->edge->from);
      vertex.insert(p->edge->to);
    }
    std::reverse(begin(), end());
  }
}


std::string History::to_string() const
{
  std::stringstream ss;

  //ostream_iterator<
  for(int i = 0; i < size(); i++) {
    ss << at(i)->to_string() << "; ";
  }
  return ss.str();
}

void EmbVector::build(const Graph &graph, Emb *e)
{
  // first build history
  if(e) {
    push_back(e->edge);
    edge.insert((*e).edge.id);
    vertex.insert((*e).edge.from);
    vertex.insert((*e).edge.to);
    //cout<<e<<endl;
    for(Emb *p = e->prev; p; p = p->prev) {
      //cout<<p<<endl;
      push_back(p->edge);
      edge.insert((*p).edge.id);
      vertex.insert((*p).edge.from);
      vertex.insert((*p).edge.to);
    }
    std::reverse(begin(), end());
  }
}

void EmbVector::build(const Graph &graph, const std::vector<Embeddings2> &emb, int emb_col, int index)
{
  //build history of the embedding backwards from index of the last vector
  //std::cout<<"in history build"<<std::endl;
  for(int k = emb_col; k >= 0; k--) {
    //std::cout<< " k = "<<k <<" index = "<<index<<" emb size = "<<emb.size()<<std::endl;
    Edge e = emb[k][index].edge;
    //cout<<e.to_string() << " ";
    push_back(e);
    edge.insert(e.id);
    vertex.insert(e.from);
    vertex.insert(e.to);
    index = emb[k][index].prev;
  }
  //cout<<endl;
  std::reverse(begin(), end());
}


std::string EmbVector::to_string() const
{
  std::stringstream ss;
  //ostream_iterator<
  for(int i = 0; i < size(); i++) {
    ss << at(i).to_string() << "; ";
  }
  return ss.str();
}

std::string EmbVector::to_string_global_vid(Graph &g) const
{
  std::stringstream ss;
  //ostream_iterator<
  for(int i = 0; i < size(); i++) {
    ss << "e(" << g[(*this)[i].from].global_vid << "," << g[(*this)[i].to].global_vid << "," << (*this)[i].elabel  << ");";
  }
  return ss.str();
}

/* Original comment:
 * get_forward_pure ()
 *  e1 (from1, elabel1, to1)
 *  from edge e2(from2, elabel2, to2)
 *
 *  minlabel <= elabel2,
 *  (elabel1 < elabel2 ||
 *  (elabel == elabel2 && tolabel1 < tolabel2)
 *  (elabel1, to1)
 *
 * RK comment:
 * ???? gets the edge that starts and extends the right-most path.
 *
 */
bool get_forward_rmpath(const Graph &graph, Edge *e, int minlabel, History& history, types::EdgeList &result)
{
  result.clear();
  assert(e->to >= 0 && e->to < graph.size());
  assert(e->from >= 0 && e->from < graph.size());
  //if(e->to >= graph.size())
  //cout<< " e->from " << e->from <<  " e->to " << e->to <<endl;
  //if(e->from >= graph.size())
  //cout<< " e->from " << e->from <<  " e->to " << e->to <<endl;

  int tolabel = graph[e->to].label;

  for(Vertex::const_edge_iterator it = graph[e->from].edge.begin();
      it != graph[e->from].edge.end(); ++it) {
    int tolabel2 = graph[it->to].label;
    //if(it->to >= graph.size())
    //cout<< " it->to " << it->to <<endl;

    if(e->to == it->to || minlabel > tolabel2 || history.hasVertex(it->to))
      continue;

    if(e->elabel < it->elabel || (e->elabel == it->elabel && tolabel <= tolabel2))
      result.push_back(const_cast<Edge*>(&(*it)));
  }

  return (!result.empty());
}

bool get_forward_rmpath(const Graph &graph, Edge e, int minlabel, EmbVector& history, std::vector<Edge> &result)
{
  result.clear();
  assert(e.to >= 0 && e.to < graph.size());
  assert(e.from >= 0 && e.from < graph.size());
  //if(e.to >= graph.size())
  //cout<< " e.from " << e.from <<  " e.to " << e.to <<endl;
  //if(e.from >= graph.size())
  //cout<< " e.from " << e.from <<  " e.to " << e.to <<endl;

  int tolabel = graph[e.to].label;

  for(Vertex::const_edge_iterator it = graph[e.from].edge.begin();
      it != graph[e.from].edge.end(); ++it) {
    int tolabel2 = graph[it->to].label;
    //if(it->to >= graph.size())
    //cout<< " it->to " << it->to <<endl;

    if(e.to == it->to || minlabel > tolabel2 || history.hasVertex(it->to))
      continue;

    if(e.elabel < it->elabel || (e.elabel == it->elabel && tolabel <= tolabel2))
      result.push_back(*it);
  }

  return (!result.empty());
}


/* Original comment:
 * get_forward_pure ()
 *  e (from, elabel, to)
 * RK comment: this function takes a "pure" forward edge, that is: an
 * edge that extends the last node of the right-most path, i.e., the
 * right-most node.
 *
 */
bool get_forward_pure(const Graph &graph, Edge *e, int minlabel, History& history, types::EdgeList &result)
{
  result.clear();

  assert(e->to >= 0 && e->to < graph.size());
  //if(e->to >= graph.size())
  //cout<< " e->from " << e->from <<  " e->to " << e->to <<endl;

  // Walk all edges leaving from vertex e->to.
  for(Vertex::const_edge_iterator it = graph[e->to].edge.begin();
      it != graph[e->to].edge.end(); ++it) {
    // -e-> [e->to] -it-> [it->to]
    assert(it->to >= 0 && it->to < graph.size());
    //if(it->to >= graph.size())
    //cout<< " it->to " << it->to <<endl;

    if(minlabel > graph[it->to].label || history.hasVertex(it->to))
      continue;

    result.push_back(const_cast<Edge*>(&(*it)));
  }

  return (!result.empty());
}

bool get_forward_pure(const Graph &graph, Edge e, int minlabel, EmbVector& history, std::vector<Edge> &result)
{
  result.clear();

  assert(e.to >= 0 && e.to < graph.size());
  //if(e.to >= graph.size())
  //cout<< " e.from " << e.from <<  " e.to " << e.to <<endl;

  // Walk all edges leaving from vertex e->to.
  for(Vertex::const_edge_iterator it = graph[e.to].edge.begin();
      it != graph[e.to].edge.end(); ++it) {
    // -e-> [e->to] -it-> [it->to]
    assert(it->to >= 0 && it->to < graph.size());
    //if(it->to >= graph.size())
    //cout<< " it->to " << it->to <<endl;

    if(minlabel > graph[it->to].label || history.hasVertex(it->to))
      continue;

    result.push_back(*it);
  }

  return (!result.empty());
}



bool get_forward_root(const Graph &g, const Vertex &v, types::EdgeList &result)
{
  result.clear();
  for(Vertex::const_edge_iterator it = v.edge.begin(); it != v.edge.end(); ++it) {
    assert(it->to >= 0 && it->to < g.size());
    if(v.label <= g[it->to].label)
      result.push_back(const_cast<Edge*>(&(*it)));
  }

  return (!result.empty());
}

bool get_forward_root(const Graph &g, const Vertex &v, std::vector<Edge> &result)
{
  result.clear();
  for(Vertex::const_edge_iterator it = v.edge.begin(); it != v.edge.end(); ++it) {
    assert(it->to >= 0 && it->to < g.size());
    if(v.label <= g[it->to].label)
      result.push_back(*it);
  }
  return (!result.empty());
}


/* Original comment:
 *  get_backward (graph, e1, e2, history);
 *  e1 (from1, elabel1, to1)
 *  e2 (from2, elabel2, to2)
 *  to2 -> from1
 *
 *  (elabel1 < elabel2 ||
 *  (elabel == elabel2 && tolabel1 < tolabel2) . (elabel1, to1)
 *
 * RK comment: gets backward edge that starts and ends at the right most path
 * e1 is the forward edge and the backward edge goes to e1->from
 */
Edge *get_backward(const Graph &graph, Edge* e1, Edge* e2, History& history)
{
  if(e1 == e2)
    return 0;

  assert(e1->from >= 0 && e1->from < graph.size());
  assert(e1->to >= 0 && e1->to < graph.size());
  assert(e2->to >= 0 && e2->to < graph.size());

  for(Vertex::const_edge_iterator it = graph[e2->to].edge.begin();
      it != graph[e2->to].edge.end(); ++it) {
    if(history.hasEdge(it->id))
      continue;

    if((it->to == e1->from) &&
       ((e1->elabel < it->elabel) ||
        (e1->elabel == it->elabel) &&
        (graph[e1->to].label <= graph[e2->to].label)
       )) {
      return const_cast<Edge*>(&(*it));
    } // if(...)
  } // for(it)

  return 0;
}

bool get_backward(const Graph &graph, Edge e1, Edge e2, EmbVector& history, Edge& result)
{
  if(e1.from == e2.from && e1.to == e2.to && e1.elabel == e2.elabel)
    return false;

  assert(e1.from >= 0 && e1.from < graph.size());
  assert(e1.to >= 0 && e1.to < graph.size());
  assert(e2.to >= 0 && e2.to < graph.size());

  for(Vertex::const_edge_iterator it = graph[e2.to].edge.begin();
      it != graph[e2.to].edge.end(); ++it) {
    if(history.hasEdge(*it))
      continue;

    if((it->to == e1.from) &&
       ((e1.elabel < it->elabel) ||
        (e1.elabel == it->elabel) &&
        (graph[e1.to].label <= graph[e2.to].label)
       )) {

      result = *it;
      return true;
    } // if(...)
  } // for(it)

  return false;
}




std::string Graph::to_string() const
{
  std::stringstream ss;
  for(int i = 0; i < vertex_size(); i++) {
    const Vertex &v = at(i);
    for(int k = 0; k < v.edge.size(); k++) {
      //if (v.edge[k].from < v.edge[k].to){
      //ss << "from: " << v.edge[k].from << "; to: " << v.edge[k].to
      ss << "from: " << v.edge[k].from << " ( " << (*this)[v.edge[k].from].global_vid << ") "
      << "; to: " << v.edge[k].to << " ( " << (*this)[v.edge[k].to].global_vid << ") "
      << "; (" << v.label << ", " << v.edge[k].elabel << ", " << get_vertex_label(v.edge[k].to) << ")" << std::endl;
      //<< get_vertex_label(v.edge[k].from) << ", " << get_vertex_label(v.edge[k].to)<< std::endl;
      //}
    } // for k
  } // for i
  return ss.str();

  //DFSCode dfs = get_min_dfs_code();
  //return dfs.to_string();
}

////////////////////////////////////////////////////////////////////////////
//
// Obtain forward and backward extensions in the graph, given the DFS code
// (DFS_CODE and history have one to one mapping of the edges)
//
///////////////////////////////////////////////////////////////////////////

bool  get_forward(const Graph &graph, const types::DFSCode &DFS_CODE, History& history, types::EdgeList &result){

  result.clear();

  //forward extenstion from dfs_from <=> from
  int dfs_from = DFS_CODE.back().from;
  int from;

  //skip the last one in dfs code
  // get the "from" vertex id from the history
  for(int i = DFS_CODE.size() - 2; i >= 0; i-- ) {
    if( dfs_from == DFS_CODE[i].from) {
      from = history[i]->from;
      break;
    }

    if( dfs_from == DFS_CODE[i].to) {
      from = history[i]->to;
      break;
    }

  }

  types::DFS dfs = DFS_CODE.back();

  for(Vertex::const_edge_iterator it = graph[from].edge.begin(); it != graph[from].edge.end(); ++it) {
    if( it->elabel == dfs.elabel && graph[it->to].label == dfs.tolabel &&  !history.hasVertex(it->to) )
      result.push_back(const_cast<Edge*>(&(*it)));
  }

  return (!result.empty());
}

bool get_forward(const Graph &graph, const types::DFSCode &DFS_CODE, EmbVector& history, std::vector<Edge> &result){

  result.clear();

  //forward extenstion from dfs_from <=> from
  int dfs_from = DFS_CODE.back().from;
  int from;

  //skip the last one in dfs code
  // get the "from" vertex id from the history
  for(int i = DFS_CODE.size() - 2; i >= 0; i-- ) {
    if( dfs_from == DFS_CODE[i].from) {
      from = history[i].from;
      break;
    }

    if( dfs_from == DFS_CODE[i].to) {
      from = history[i].to;
      break;
    }

  }

  types::DFS dfs = DFS_CODE.back();

  for(Vertex::const_edge_iterator it = graph[from].edge.begin(); it != graph[from].edge.end(); ++it) {
    if( it->elabel == dfs.elabel && graph[it->to].label == dfs.tolabel &&  !history.hasVertex(it->to) )
      result.push_back(*it);
  }

  return (!result.empty());
}



Edge *get_backward(const Graph &graph, const types::DFSCode &DFS_CODE,  History& history)
{

  std::map<int, int> vertex_id_map;
  for(int i = 0; i<history.size(); i++) {
    if(vertex_id_map.count(DFS_CODE[i].from) == 0)
      vertex_id_map[DFS_CODE[i].from] = history[i]->from;
    if(vertex_id_map.count(DFS_CODE[i].to) == 0)
      vertex_id_map[DFS_CODE[i].to]   = history[i]->to;
  }

  //now add the backward edge using the last entry of the DFS code
  int from = vertex_id_map[DFS_CODE.back().from];
  int to   = vertex_id_map[DFS_CODE.back().to];

  for(Vertex::const_edge_iterator it = graph[from].edge.begin(); it != graph[from].edge.end(); ++it) {

    if(it->to == to)
      return const_cast<Edge*>(&(*it));

  } // for(it)

  return 0;
}

bool get_backward(const Graph &graph, const types::DFSCode &DFS_CODE, EmbVector& history, Edge& e)
{
  std::map<int, int> vertex_id_map;
  for(int i = 0; i<history.size(); i++) {
    if(vertex_id_map.count(DFS_CODE[i].from) == 0)
      vertex_id_map[DFS_CODE[i].from] = history[i].from;
    if(vertex_id_map.count(DFS_CODE[i].to) == 0)
      vertex_id_map[DFS_CODE[i].to]   = history[i].to;
  }

  //now add the backward edge using the last entry of the DFS code
  int from = vertex_id_map[DFS_CODE.back().from];
  int to   = vertex_id_map[DFS_CODE.back().to];

  for(Vertex::const_edge_iterator it = graph[from].edge.begin(); it != graph[from].edge.end(); ++it) {

    if(it->to == to) {
      e = *it;
      return true;
    }

  } // for(it)

  return false;
}

void Graph::delete_edge(int from, int to){

  for(Vertex::edge_iterator it = (*this)[from].edge.begin(); it != (*this)[from].edge.end(); ) {
    if(it->to == to) {
      (*this)[from].edge.erase(it);
    }else{
      ++it;
    }
  }

  for(Vertex::edge_iterator it = (*this)[to].edge.begin(); it != (*this)[to].edge.end(); ) {
    if(it->to == from) {
      (*this)[to].edge.erase(it);
    }else{
      ++it;
    }
  }

  std::vector<int> d;
  if((*this)[from].edge.size() == 0) d.push_back(from);
  if((*this)[to].edge.size() == 0) d.push_back(to);

  if(d.size() > 0)
    delete_vertices(d);
}

//only applicable for adjp/ladjp format
void Graph::delete_vertices(std::vector<int> local_ids){

  if(local_ids.size() == 0) return;

  std::set<int> to_delete_ids_set;
  for(std::vector<int>::iterator it = local_ids.begin(); it != local_ids.end(); it++) {
    to_delete_ids_set.insert(*it);
  }

  //find the vertices that are only connected to "to delete" vertices, put them in local_ids list
  for(int vid = 0; vid < this->size(); vid++) {

    if(to_delete_ids_set.count(vid) > 0)
      continue;

    bool all_neighbors_deleted = true;
    for(std::vector<Edge>::iterator edge_iterator = (*this)[vid].edge.begin();
        edge_iterator != (*this)[vid].edge.end(); edge_iterator++ ) {
      if(to_delete_ids_set.count(edge_iterator->to) == 0) {
        all_neighbors_deleted = false;
        break;
      }
    }
    if(all_neighbors_deleted) {
      local_ids.push_back(vid);
      if(vid <= max_local_vid)
        max_local_vid--;
    }

  }


  //create old -> new local id map, excluding the deleted vertices (local_ids)
  //also update global_local_id map
  std::sort(local_ids.begin(), local_ids.end());
  std::map<int, int> new_local_ids;
  int i = 0, new_id = 0;
  for(int lid = 0; lid< this->size(); lid++) {
    // if the vertex is in delete list remove from global local id map
    if(local_ids[i] == lid) {
      global_local_id_map.erase((*this)[lid].global_vid);
      //skip and advance to the next vertex
      i++;
    }else{
      //else store the new id and update global local id map
      new_local_ids[lid] = new_id;                   //old_id -> new_id
      global_local_id_map[(*this)[lid].global_vid] = new_id;
      new_id++;
    }
  }

  //cout<<utils::print_map(global_local_id_map)<<endl;

  //update the neighbor lists, remove unnecessary neighbors, update neighbors with their new local ids
  for(int vid = 0; vid < this->size(); vid++) {

    if(new_local_ids.count(vid) > 0 ) {
      //for(int i=0; i < (*this)[vid].edge.size(); i++){
      for(std::vector<Edge>::iterator edge_iterator = (*this)[vid].edge.begin();
          edge_iterator != (*this)[vid].edge.end(); ) {
        //int to = (*this)[vid].edge[i].to;
        edge_iterator->from = new_local_ids[vid];
        if(new_local_ids.count(edge_iterator->to) == 0) {
          //remove it
          (*this)[vid].edge.erase(edge_iterator);
        }else{
          //update it
          edge_iterator->to = new_local_ids[edge_iterator->to];
          ++edge_iterator;
        }
      }



    } else {
      //to be deleted vertices, no need to update, will be deleted in next step
    }

  }

  //Finally erase the vertices from the graph
  int vid = 0;
  for(std::vector<Vertex>::iterator it = this->begin() + vid; it != this->end(); ) {
    if(new_local_ids.count(vid) == 0 ) {
      this->erase(it);
    }else{
      ++it;
    }
    vid++;
  }

  //cout<<" number of vertices = "<<this->size()<<endl;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Serialization
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


size_t Vertex::get_serialized_size(const Vertex &vrtx)
{
  //    vertex label + 4 * #of edges  * sizeof(int)    +   number of edges  + label;
  return sizeof(int) + 4 * vrtx.edge.size() * sizeof(int) + sizeof(int)     + sizeof(int);
}


size_t Vertex::get_serialized_size(char *buffer, size_t buffer_size)
{
  int s = *((int*)buffer);
  return s;
}


size_t Vertex::serialize(const Vertex &vrtx, char *buffer, size_t buffer_size)
{
  if(buffer_size < get_serialized_size(vrtx)) throw std::runtime_error("Buffer too small.");
  int pos = 0;

  // size of this serialized vertex in bytes.
  *((int*)(buffer + pos)) = get_serialized_size(vrtx);
  pos += sizeof(int);

  // store the vertex label
  *((int*)(buffer + pos)) = vrtx.label;
  pos += sizeof(int);


  // store number of edges
  *((int*)(buffer + pos)) = vrtx.edge.size();
  pos += sizeof(int);

  for(int i = 0; i < vrtx.edge.size(); i++) {
    *((int*)(buffer + pos)) = vrtx.edge[i].from;
    pos += sizeof(int);

    *((int*)(buffer + pos)) = vrtx.edge[i].to;
    pos += sizeof(int);

    *((int*)(buffer + pos)) = vrtx.edge[i].elabel;
    pos += sizeof(int);

    *((int*)(buffer + pos)) = vrtx.edge[i].id;
    pos += sizeof(int);
  } // for i

  return pos;
} // Vertex::serialize


size_t Vertex::deserialize(Vertex &vrtx, char *buffer, size_t buffer_size)
{
  // TODO: check minimum buffer size
  if(buffer_size < get_serialized_size(buffer, buffer_size)) throw std::runtime_error("Buffer too small.");
  int pos = 0;
  vrtx.edge.clear();

  // read buffer s
  pos += sizeof(int);

  // read the vertex label
  vrtx.label = *((int*)(buffer + pos));
  pos += sizeof(int);

  // read the number of edges
  int edge_count = *((int*)(buffer + pos));
  pos += sizeof(int);


  for(int i = 0; i < edge_count; i++) {
    Edge tmp_edge;
    tmp_edge.from = *((int*)(buffer + pos));
    pos += sizeof(int);
    tmp_edge.to = *((int*)(buffer + pos));
    pos += sizeof(int);
    tmp_edge.elabel = *((int*)(buffer + pos));
    pos += sizeof(int);
    tmp_edge.id = *((int*)(buffer + pos));
    pos += sizeof(int);
    vrtx.edge.push_back(tmp_edge);
  } // for i

  return pos;
} // Vertex::deserialize

size_t Graph::get_serialized_size_for_partition(int vid)
{
  //    #of edges + 4 * #of edges  * sizeof(int) ;
  return sizeof(int) + 4 * (*this)[vid].edge.size() * sizeof(int);
}

size_t Graph::serialize_neighbors_for_partition(int global_vid, int *buffer, int buffer_size)
{
  int local_vid = (*this).global_local_id_map[global_vid];

  if(local_vid < 0) throw std::runtime_error("invalid global vid");

  if(buffer_size < get_serialized_size_for_partition(local_vid)) throw std::runtime_error("Buffer too small.");
  int pos = 0;
  int num_edges = 0;

  // store number of edges not pertaining to the requesting partition

  *((int*)(buffer + pos)) = (*this)[local_vid].edge.size();
  pos += sizeof(int);


  for(int i = 0; i < (*this)[local_vid].edge.size(); i++) {

    *((int*)(buffer + pos)) = (*this)[local_vid].edge[i].elabel;
    pos += sizeof(int);

    *((int*)(buffer + pos)) = (*this)[(*this)[local_vid].edge[i].to].global_vid;
    pos += sizeof(int);

    *((int*)(buffer + pos)) = (*this)[(*this)[local_vid].edge[i].to].label;
    pos += sizeof(int);

    *((int*)(buffer + pos)) = (*this)[(*this)[local_vid].edge[i].to].orig_part_id;
    pos += sizeof(int);
  } // for i

  return pos;
} // Graph::serialize


size_t Graph::get_serialized_size_for_partition(int vid, int requester_partition_id)
{
  int num_edges = 0;
  for(int i = 0; i < (*this)[vid].edge.size(); i++)
    if((*this)[(*this)[vid].edge[i].to].orig_part_id != requester_partition_id )
      num_edges++;
  // from node + #of edges + 4 * #of edges * sizeof(int) ;
  //cout <<  "serialized size =" << (sizeof(int) + sizeof(int) + 4 * num_edges * sizeof(int)) <<endl;
  return sizeof(int) + sizeof(int) + 4 * num_edges * sizeof(int);
}

size_t Graph::serialize_neighbors_for_partition(int requester_partition_id, int global_vid, int* &buffer, int buffer_size)
{
  int local_vid = (*this).global_local_id_map[global_vid];
  if(local_vid < 0) throw std::runtime_error("invalid global vid");

  //if(buffer_size < get_serialized_size_for_partition(local_vid)) throw std::runtime_error("Buffer too small.");
  int pos = 0;
  int num_edges = 0;

  // store number of edges not pertaining to the requesting partition
  for(int i = 0; i < (*this)[local_vid].edge.size(); i++)
    if((*this)[(*this)[local_vid].edge[i].to].orig_part_id != requester_partition_id )
      num_edges++;


  //*((int*)(buffer + pos)) = global_vid;
  //pos += sizeof(int);
  buffer[pos++] = global_vid;

  buffer[pos++] = num_edges;

  for(int i = 0; i < (*this)[local_vid].edge.size(); i++) {
    int to_local_vid = (*this)[local_vid].edge[i].to;
    if((*this)[to_local_vid].orig_part_id != requester_partition_id ) {

      //if( (*this)[to_local_vid].global_vid == -1)
      //cout<<" to local id = " <<to_local_vid << " global id "<<(*this)[to_local_vid].global_vid <<" graph size= "<< this->size() <<endl;
      buffer[pos++] = (*this)[to_local_vid].global_vid;
      buffer[pos++] = (*this)[to_local_vid].label;
      buffer[pos++] = (*this)[to_local_vid].orig_part_id;
      buffer[pos++] = (*this)[local_vid].edge[i].elabel;
    }
  } // for i

  return pos;
} // Graph::serialize

size_t Graph::serialize_neighbors_for_partition(int requester_partition_id, int global_vid, int* &buffer, int buffer_size, std::set<int> &exclusions)
{
  int local_vid = (*this).global_local_id_map[global_vid];
  if(local_vid < 0) throw std::runtime_error("invalid global vid");

  //if(buffer_size < get_serialized_size_for_partition(local_vid)) throw std::runtime_error("Buffer too small.");
  int pos = 0;
  int num_edges = 0;

  // store number of edges not pertaining to the requesting partition
  for(int i = 0; i < (*this)[local_vid].edge.size(); i++)
    if((*this)[(*this)[local_vid].edge[i].to].orig_part_id != requester_partition_id
       && exclusions.find((*this)[(*this)[local_vid].edge[i].to].global_vid) == exclusions.end() )
      num_edges++;


  //*((int*)(buffer + pos)) = global_vid;
  //pos += sizeof(int);
  buffer[pos++] = global_vid;

  buffer[pos++] = num_edges;

  for(int i = 0; i < (*this)[local_vid].edge.size(); i++) {
    int to_local_vid = (*this)[local_vid].edge[i].to;
    if((*this)[to_local_vid].orig_part_id != requester_partition_id
       && exclusions.find((*this)[(*this)[local_vid].edge[i].to].global_vid) == exclusions.end() ) {

      //if( (*this)[to_local_vid].global_vid == -1)
      //cout<<" to local id = " <<to_local_vid << " global id "<<(*this)[to_local_vid].global_vid <<" graph size= "<< this->size() <<endl;
      buffer[pos++] = (*this)[to_local_vid].global_vid;
      buffer[pos++] = (*this)[to_local_vid].label;
      buffer[pos++] = (*this)[to_local_vid].orig_part_id;
      buffer[pos++] = (*this)[local_vid].edge[i].elabel;
    }
  } // for i

  return pos;
} // Graph::serialize



size_t Graph::deserialize_neighbors_for_partition(int part_id, int*&buffer, int buffer_size, int &global_vid_recv)
{
  int pos = 0;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // read the global id
  int global_vid = buffer[pos++];

  global_vid_recv = global_vid;
  int from_local_vid = this->global_local_id_map[global_vid];
  //make the from vertex a local vertex to the same partition, since it now has all the neighbors
  (*this)[from_local_vid].vertex_part_id = part_id;

  // read the number of edges
  int edge_count = buffer[pos++];

  for(int i = 0; i < edge_count; i++) {
    int to_global_vid = buffer[pos++];
    int to_local_vid;
    if( global_local_id_map.count(to_global_vid) == 0) {
      to_local_vid = this->size();
      this->resize(to_local_vid + 1);
      //if(part_id == 0)
      //	std::cout << "rank " << rank << " === SIZE CHANGED from " << to_local_vid << " to " <<this->size()<<std::endl;
      global_local_id_map[to_global_vid] = to_local_vid;
      (*this)[to_local_vid].global_vid = to_global_vid;
      (*this)[to_local_vid].label = buffer[pos++];
      (*this)[to_local_vid].vertex_part_id = buffer[pos++];
      (*this)[to_local_vid].orig_part_id = (*this)[to_local_vid].vertex_part_id;
    } else {
      to_local_vid = global_local_id_map[to_global_vid];
      pos += 2;           // label and part_id are already there
    }
    //if(part_id == 0)
    //std::cout << "rank " << rank << " Graph modified: from_local_vid  "<< from_local_vid << " to_local_vid " << to_local_vid << " label " << (*this)[to_local_vid].label << " vertex_part_id " << (*this)[to_local_vid].vertex_part_id <<endl;
    int elabel = buffer[pos++];

    types::Edge r;
    //eliminate duplicate here...
    if(!(*this)[from_local_vid].find(from_local_vid, to_local_vid, r)) {
      (*this)[from_local_vid].push(from_local_vid, to_local_vid, elabel);
      if(!directed)
        (*this)[to_local_vid].push(to_local_vid, from_local_vid, elabel);
    }
    //this->buildEdge();

  } // for i

  return pos;
} // Vertex::deserialize


size_t Graph::deserialize_multiple_neighbors_for_partition(int part_id, int*&buffer, int buffer_size, int num_neighbors)
{
  int pos = 0;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //std::stringstream ss;
  //ss << "rank "<<rank << " received neighbor data ";
  for(int k = 0; k < num_neighbors; k++) {
    // read the global id
    int global_vid = buffer[pos++];
    //ss << global_vid << " ";
    int from_local_vid = this->global_local_id_map[global_vid];
    //make the from vertex a local vertex to the same partition, since it now has all the neighbors
    (*this)[from_local_vid].vertex_part_id = part_id;

    // read the number of edges
    int edge_count = buffer[pos++];

    for(int i = 0; i < edge_count; i++) {
      int to_global_vid = buffer[pos++];
      int to_local_vid;
      if( global_local_id_map.count(to_global_vid) == 0) {
        to_local_vid = this->size();
        this->resize(to_local_vid + 1);
        //if(part_id == 0)
        //	std::cout << "rank " << rank << " === SIZE CHANGED from " << to_local_vid << " to " <<this->size()<<std::endl;
        global_local_id_map[to_global_vid] = to_local_vid;
        (*this)[to_local_vid].global_vid = to_global_vid;
        (*this)[to_local_vid].label = buffer[pos++];
        (*this)[to_local_vid].vertex_part_id = buffer[pos++];
        (*this)[to_local_vid].orig_part_id = (*this)[to_local_vid].vertex_part_id;
      } else {
        to_local_vid = global_local_id_map[to_global_vid];
        pos += 2;                 // label and part_id are already there
      }
      //if(part_id == 0)
      //std::cout << "rank " << rank << " Graph modified: from_local_vid  "<< from_local_vid << " to_local_vid " << to_local_vid << " label " << (*this)[to_local_vid].label << " vertex_part_id " << (*this)[to_local_vid].vertex_part_id <<endl;
      int elabel = buffer[pos++];

      types::Edge r;
      //eliminate duplicate here...
      if(!(*this)[from_local_vid].find(from_local_vid, to_local_vid, r)) {
        (*this)[from_local_vid].push(from_local_vid, to_local_vid, elabel);
        if(!directed)
          (*this)[to_local_vid].push(to_local_vid, from_local_vid, elabel);
      }
    }       // for i
  }
  //if(num_neighbors >0)
  //std::cout<<ss.str()<<endl;
  return pos;

}


size_t Graph::get_serialized_size(const Graph &grph)
{
  size_t s = sizeof(int) + sizeof(int) + sizeof(int) + sizeof(bool); // edge_size_ + total buffer size + number of vertices + variable directed(bool)
  for(int i = 0; i < grph.size(); i++) {
    s += Vertex::get_serialized_size(grph[i]);
  } // for i
  return s;
} // Graph::get_serialized_size


size_t Graph::get_serialized_size(char *buffer, size_t buffer_size)
{
  return *((int*) buffer);
}


size_t Graph::serialize(const Graph &grph, char *buffer, size_t buffer_size)
{
  if(get_serialized_size(grph) > buffer_size) throw std::runtime_error("Buffer too small.");
  int pos = 0;

  // store buffer size
  *((int*)(buffer + pos)) = get_serialized_size(grph);
  pos += sizeof(int);

  // store edge_size_
  *((int*)(buffer + pos)) = grph.edge_size_;
  pos += sizeof(int);

  // store number of vertices
  *((bool*)(buffer + pos)) = grph.directed;
  pos += sizeof(grph.directed);


  // store number of vertices
  *((int*)(buffer + pos)) = grph.size();
  pos += sizeof(int);

  for(int i = 0; i < grph.size(); i++) {
    int tmp_pos = Vertex::serialize(grph.at(i), buffer + pos, buffer_size - pos);
    pos += tmp_pos;
  } // for i

  return pos;
} // Graph::serialize


size_t Graph::deserialize(Graph &grph, char *buffer, size_t buffer_size)
{
  if(Graph::get_serialized_size(buffer, buffer_size) > buffer_size) throw std::runtime_error("Buffer too small.");

  grph.clear();
  int pos = 0;

  // store buffer size
  pos += sizeof(int);

  // store edge_size_
  grph.edge_size_ = *((int*)(buffer + pos));
  pos += sizeof(int);

  // store number of vertices
  grph.directed = *((bool*)(buffer + pos));
  pos += sizeof(grph.directed);


  // store number of vertices
  int vert_count = *((int*)(buffer + pos));
  pos += sizeof(int);

  for(int i = 0; i < vert_count; i++) {
    Vertex tmp_vert;
    int tmp_pos = Vertex::deserialize(tmp_vert, buffer + pos, buffer_size - pos);
    grph.push_back(tmp_vert);
    pos += tmp_pos;
  } // for i

  return pos;
} // Graph::deserialize





size_t Graph::get_serialized_size(const graph_database_t &grph_db)
{
  size_t min_buff_size = 0;

  min_buff_size += sizeof(int) + sizeof(int); // size of the database + size of the buffer

  for(size_t i = 0; i < grph_db.size(); i++) {
    min_buff_size += get_serialized_size(grph_db[i]);
  } // for i

  return min_buff_size;
} // Graph::get_serialized_size


size_t Graph::get_serialized_size_db(char *buffer, size_t buffer_size)
{
  //abort();
  return *((int*) buffer);
} // Graph::get_serialized_size


size_t Graph::serialize(const graph_database_t &grph_db, char *buffer, size_t buffer_size)
{
  size_t pos = 0;

  int min_buff_size = get_serialized_size(grph_db);
  if(min_buff_size > buffer_size) throw std::runtime_error("Buffer too small.");

  *((int*)(buffer + pos)) = min_buff_size;
  pos += sizeof(int);

  *((int*)(buffer + pos)) = grph_db.size();
  pos += sizeof(int);

  for(int i = 0; i < grph_db.size(); i++) {
    size_t tmp_pos = serialize(grph_db[i], buffer + pos, buffer_size - pos);
    pos += tmp_pos;
  }

  return pos;
} // Graph::serialize


size_t Graph::deserialize(graph_database_t &grph_db, char *buffer, size_t buffer_size)
{
  int min_buf_size = get_serialized_size_db(buffer, buffer_size);

  if(buffer_size < min_buf_size) throw std::runtime_error("Buffer too small.");

  grph_db.clear();

  size_t pos = 0;
  // skip buffer size
  pos += sizeof(int);

  int grph_db_size = *((int*)(buffer + pos));
  pos += sizeof(int);

  for(int i = 0; i < grph_db_size; i++) {
    Graph grph;
    size_t tmp_pos = deserialize(grph, buffer + pos, buffer_size - pos);
    pos += tmp_pos;
    grph_db.push_back(grph);
  } // for i

  return pos;
} // Graph::deserialize


} // namespace types


