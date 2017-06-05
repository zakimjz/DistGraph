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



#include <external_neighbor_handler.hpp>
#include <unistd.h>
#include <alloc_tools.hpp>

using namespace std;

namespace algs {

void external_neighbor_handler::init_external_neighbor_handler(unsigned int rank, unsigned int numtasks, unsigned int numthreads)
{
  this->rank = rank;
  this->numtasks = numtasks;
  this->numthreads = numthreads;

  this->num_partitions = numtasks;



  for(int i = 0; i<numthreads; i++) {
    std::map<std::string, std::set<int> > t;
    threads_waiting.push_back(false);
    requests_for_dfscodes_threads.push_back(t);
  }

  for(int i = 0; i < numtasks; i++) {
    std::set<int> t;
    processed_external_neighbors[i] = t;
  }

  computation_finished = false;
  if(rank == 0)
    has_token = true;
  else
    has_token = false;

  logger = Logger::get_logger("EXTERNAL_NEIGHBOR_HANDLER");
  TRACE(*logger, "external neighbor handler initialized");
}

void external_neighbor_handler::set_num_partitions(int num_part){
  num_partitions = num_part;
}

/*
 * returns: -1 if no request is made, otherwise rank of the requesting processor
 */
int external_neighbor_handler::check_message(int &size)
{
  MPI_Status stat;
  int flag;
  TRACE5(*logger, "checking message");
  MPI_Iprobe(MPI_ANY_SOURCE, MSG_NEIGH_DATA, MPI_COMM_WORLD, &flag, &stat);
  DEBUG(*logger, "prob value =" << flag );

  if(flag) {
    MPI_Get_count(&stat, MPI_INT, &size);
    DEBUG(*logger, "got request from: " << stat.MPI_SOURCE);
    return stat.MPI_SOURCE;
  }
  return -1;
}

void external_neighbor_handler::union_neighbor_requests_from_rank_group(){
  MPI_Comm rank_group;
  int membershipKey = rank % num_partitions;
  MPI_Comm_split(MPI_COMM_WORLD, membershipKey, rank, &rank_group);

  int group_size = numtasks / num_partitions;
  int send_length = 2 * external_neighbor_requests.size();
  int *recv_length = new int[group_size];

  all_gather_msg(&send_length, recv_length, 1, rank_group);

  DEBUG(*logger, "received lengths = " << utils::print_array(recv_length, group_size));

  int total_recv_length = 0;
  for(int i = 0; i < group_size; i++)
    total_recv_length += recv_length[i];

  int *send_buf, *recv_buf;
  //int *send_buf = new int[send_length];
  //int *recv_buf = new int[total_recv_length];
  NEW_INT_ARRAY(send_buf, send_length, *logger);
  NEW_INT_ARRAY(recv_buf, total_recv_length, *logger);

  int offset = 0;
  for(std::map<int, int>::iterator it = external_neighbor_requests.begin(); it != external_neighbor_requests.end(); ++it) {
    send_buf[offset++] = it->first;
    send_buf[offset++] = it->second;
  }

  all_gatherv_msg(send_buf, recv_buf, send_length, recv_length, rank_group);

  DEBUG(*logger, utils::print_array(recv_buf, total_recv_length));

  offset = 0;
  for(int i = 0; i < group_size; i++) {
    if(i == rank / num_partitions ) {
      offset += recv_length[i];
    } else {
      for(int j = 0; j < recv_length[i] / 2; j++) {
        int global_vid = recv_buf[offset++];
        int dest = recv_buf[offset++];
        if(external_neighbor_requests.count(global_vid) == 0)
          external_neighbor_requests[global_vid] = dest;
      }
    }
  }

  DEBUG(*logger, "done inserting total requests = " << external_neighbor_requests.size() );

  delete [] recv_length;
  //delete [] send_buf;
  //delete [] recv_buf;
  DELETE_INT_ARRAY(send_buf, *logger);
  DELETE_INT_ARRAY(recv_buf, *logger);

}


void external_neighbor_handler::union_used_pseudo_locals_in_rank_group(std::set<int> &used_pseudo_locals){
  MPI_Comm rank_group;
  int membershipKey = rank % num_partitions;
  MPI_Comm_split(MPI_COMM_WORLD, membershipKey, rank, &rank_group);

  int group_size = numtasks / num_partitions;
  int send_length = used_pseudo_locals.size();
  int *recv_length = new int[group_size];

  all_gather_msg(&send_length, recv_length, 1, rank_group);

  DEBUG(*logger, "received lengths = " << utils::print_array(recv_length, group_size));

  int total_recv_length = 0;
  for(int i = 0; i < group_size; i++)
    total_recv_length += recv_length[i];

  int *send_buf, *recv_buf;
  //int *send_buf = new int[send_length];
  //int *recv_buf = new int[total_recv_length];
  NEW_INT_ARRAY(send_buf, send_length, *logger);
  NEW_INT_ARRAY(recv_buf, total_recv_length, *logger);

  int offset = 0;
  for(std::set<int>::iterator it = used_pseudo_locals.begin(); it != used_pseudo_locals.end(); ++it) {
    send_buf[offset++] = *it;
  }

  all_gatherv_msg(send_buf, recv_buf, send_length, recv_length, rank_group);

  DEBUG(*logger, utils::print_array(recv_buf, total_recv_length));

  offset = 0;
  for(int i = 0; i < group_size; i++) {
    if(i == rank / num_partitions ) {
      offset += recv_length[i];
    } else{
      for(int j = 0; j < recv_length[i]; j++) {
        //if(used_pseudo_locals.count(recv_buf[offset]) == 0)
        used_pseudo_locals.insert(recv_buf[offset++]);
      }
    }
  }

  INFO(*logger, "done inserting used pseudo locals = " << used_pseudo_locals.size() );

  delete [] recv_length;
  //delete [] send_buf;
  //delete [] recv_buf;
  DELETE_INT_ARRAY(send_buf, *logger);
  DELETE_INT_ARRAY(recv_buf, *logger);
}

void external_neighbor_handler::send_recv_external_neighbors(){

  if(numtasks > num_partitions)
    union_neighbor_requests_from_rank_group();


  std::vector<std::vector<int> > requests_per_rank(num_partitions);

  for(std::map<int,int>::iterator it = external_neighbor_requests.begin(); it != external_neighbor_requests.end(); ++it) {
    int global_vid = it->first;
    int dest = it->second;
    // e.g. num_partitions = 4, numtasks = 8
    // proc 0 sends to 3 => proc 4 sends to 7 (instead of 3)
    //if(numtasks > num_partitions) {
    //dest += (rank/num_partitions) * num_partitions;
    //}
    requests_per_rank[dest].push_back(global_vid);
    DEBUG(*logger, "rank " << rank << "group id " << group_id << " sending request for vertex " << global_vid << " to " << dest);
  }

  //stringstream ss1;
  //ss1 << "rank " << rank <<" requesting number of neighbors ";
  //for(int i = 0;i < numtasks; i++)
  //ss1 << " from " << i << " : "<< requests_per_rank[i].size();

  //DEBUG(*logger, ss1.str());

  //ss1.clear();
  //ss1 << "rank " << rank <<" requesting neighbors ";
  //for(int i = 0;i < numtasks; i++)
  //ss1 << " from " << i << " : "<< utils::print_vector(requests_per_rank[i]);

  //TRACE4(*logger, ss1.str());

  //Send data only within rankgroups, for 2 partitions and 4 tasks groups will be (0,1) (2,3)
  //There will be numtasks/num_partitions groups
  MPI_Comm rank_group;
  int membershipKey = rank / num_partitions;
  MPI_Comm_split(MPI_COMM_WORLD, membershipKey, rank, &rank_group);

  int *send_lengths = new int[num_partitions];       //number of neighbor requests per rank in group
  int *recv_lengths = new int[num_partitions];

  //get neighbor requests
  for(int i = 0; i<num_partitions; i++)
    send_lengths[i] = requests_per_rank[i].size();

  alltoall_msg(send_lengths, recv_lengths, 1, rank_group);

  //get the total number of neigbor requests from all ranks
  int total_send_length = 0, total_recv_length = 0;
  for(int i = 0; i<num_partitions; i++) {
    total_send_length += send_lengths[i];
    total_recv_length += recv_lengths[i];
  }

  int *sendbuf = 0, *recvbuf = 0;
  //sendbuf = new int[total_send_length];
  //recvbuf = new int[total_recv_length];
  NEW_INT_ARRAY(sendbuf, total_send_length, *logger);
  NEW_INT_ARRAY(recvbuf, total_recv_length, *logger);

  int offset = 0;
  for(int i = 0; i<num_partitions; i++) {
    for(int j = 0; j<requests_per_rank[i].size(); j++) {
      sendbuf[offset++] = requests_per_rank[i][j];
    }
  }

  //get the global vids of the neighbor requests
  alltoallv_msg(sendbuf, recvbuf, send_lengths, recv_lengths, rank_group);

  //stringstream ss2;
  //ss2 << "rank " << rank <<" received neighbor requests ";
  offset = 0;
  for(int i = 0; i < num_partitions; i++) {
    //ss2 << " from " << i << " : "<< utils::print_array(recvbuf+offset, recv_lengths[i]);
    offset += recv_lengths[i];
  }
  //TRACE4(*logger, ss2.str());

  //Now prepare neighbor data of the received global ids
  int *send_lengths2 = new int[num_partitions];
  int *recv_lengths2 = new int[num_partitions];


  //Perform a buffered send/recv of neighbors
  //determine the send lengths
  offset = 0;
  total_send_length = 0;
  //process neighbor requests
  for(int i = 0; i<num_partitions; i++) {
    send_lengths2[i] = 1;             //num_neighbors
    recv_lengths2[i] = 0;
    for(int j = 0; j<recv_lengths[i]; j++) {
      int global_vid = recvbuf[offset++];
      int degree = get_num_neighbors(i, global_vid);
      int length = 2  + degree * 4;                   // 2 (global vid, degree) + degree*4
      send_lengths2[i] += length;
    }
    total_send_length += send_lengths2[i];
  }

  //get the buffer size for neighbor info
  alltoall_msg(send_lengths2, recv_lengths2, 1, rank_group);

  DEBUG(*logger, " send lengths  " << utils::print_array(send_lengths2, numtasks));
  total_recv_length = 0;
  for(int i = 0; i<num_partitions; i++)
    total_recv_length += recv_lengths2[i];

  //new sendbuf and recvbuf
  int *sendbuf2, *recvbuf2;
  //int *sendbuf2 = new int[total_send_length];
  //int *recvbuf2 = new int[total_recv_length];
  NEW_INT_ARRAY(sendbuf2, total_send_length, *logger);
  NEW_INT_ARRAY(recvbuf2, total_recv_length, *logger);

  //prepare the sendbufs with replies
  offset = 0;
  int send_offset = 0;
  for(int i = 0; i<num_partitions; i++) {
    sendbuf2[send_offset++] = recv_lengths[i];             //num neighbors
    for(int j = 0; j <recv_lengths[i]; j++) {
      int global_vid = recvbuf[offset++];
      send_offset += prepare_neighbor_data_to_send(i, global_vid, (sendbuf2 + send_offset), -1);
    }
  }

  stringstream ss;
  ss << " || ";
  offset = 0;
  for(int i = 0; i<num_partitions; i++) {
    ss << utils::print_array(sendbuf2 + offset, send_lengths2[i]) << " || ";
    offset += send_lengths2[i];
  }

  DEBUG(*logger, " sendbuf2 (older version) : " << ss.str());

  //send recv neighbor info
  alltoallv_msg(sendbuf2, recvbuf2, send_lengths2, recv_lengths2, rank_group);

  //process the received neighbor info
  offset = 0;
  int recv_offset = 0;
  for(int i = 0; i<num_partitions; i++) {
    receive_neighbor_data(i, recvbuf2 + recv_offset, recv_lengths2[i]);
    recv_offset += recv_lengths2[i];
  }

  ////////////////////////////


  external_neighbor_requests.clear();
  requests_for_dfscodes.clear();

  delete [] send_lengths;
  delete [] recv_lengths;

  //cout<<" sendbuf pointer " << sendbuf<<endl;

  if(sendbuf)
    DELETE_INT_ARRAY(sendbuf, *logger);
  //delete []sendbuf;
  if(recvbuf)
    DELETE_INT_ARRAY(recvbuf, *logger);
  //delete []recvbuf;

  delete [] send_lengths2;
  delete [] recv_lengths2;
  if(sendbuf2)
    DELETE_INT_ARRAY(sendbuf2, *logger);
  //delete [] sendbuf2;

  if(recvbuf2)
    DELETE_INT_ARRAY(recvbuf2, *logger);
  //delete [] recvbuf2;

}

void external_neighbor_handler::send_neighbor_requests(){

  //omp_set_lock(&lock_nr);
  int **buf = new int*[numtasks];
  std::vector<int> len;
  for(int i = 0; i < numtasks; i++) {
    //allocate maximum possible
    buf[i] = new int[1 + external_neighbor_requests.size()];
    buf[i][0] = RT_NEIGH_REQ;
    len.push_back(1);
  }

  for(std::map<int,int>::iterator it = external_neighbor_requests.begin(); it != external_neighbor_requests.end(); ++it) {
    int global_vid = it->first;
    int dest = it->second;
    // e.g. num_partitions = 4, numtasks = 8
    // proc 0 sends to 3 => proc 4 sends to 7 (instead of 3)
    if(numtasks > num_partitions) {
      dest += (rank / num_partitions) * num_partitions;
    }
    //if(is_external_neighbor_fetched.count(global_vid) == 0){
    is_external_neighbor_fetched[global_vid] = false;
    buf[dest][len[dest]++] = global_vid;
    DEBUG(*logger, "rank " << rank << "group id " << group_id << " sending request for vertex " << global_vid << " to " << dest);
    //}
  }

  for(int i = 0; i<numtasks; i++) {
    //send neighbor requests (if there is) to dest rank i
    if(len[i] > 1)
      send_msg(buf[i], len[i], i);
    delete buf[i];
  }
  delete[] buf;
  //external_neighbor_requests.clear();
  //omp_unset_lock(&lock_nr);
}

void external_neighbor_handler::send_neighbor_request(int global_vid, int dest){
  omp_set_lock(&lock_nr);
  if(is_external_neighbor_fetched.count(global_vid) == 0) {
    int buf[2];
    buf[0] = RT_NEIGH_REQ;
    buf[1] = global_vid;
    DEBUG(*logger, "rank " << rank << " sending request for vertex " << global_vid << " to " << dest);
    is_external_neighbor_fetched[global_vid] = false;
    send_msg(buf, 2, dest);
  }
  omp_unset_lock(&lock_nr);
}


bool external_neighbor_handler::receive_neighbor_data(int source, int *data, int size)
{
  if(size == 0)
    return false;
  //if(size == 3) {
  //TRACE(*logger, "Did not get any neighbors");
  //return false;
  //}
  //data++;  //skip RT_NEIGH_DATA
  MPI_Status stat;

  int num_requests = *data++;
  neighbors_fetched += num_requests;

  DEBUG(*logger, "rank " << rank << " receiving data of size=" << size << ", from " << source << " num requests " << num_requests);
  //return true;

  process_received_neighbor_data(data, size - 1, num_requests);
  data--;

  return true;
}


void external_neighbor_handler::process_neighbor_request(int source, int global_vid)
{
  int *buffer = 0;
  int degree = get_num_neighbors(source, global_vid);
  int length = 3 + degree * 4;

  DEBUG(*logger, "rank " << rank << " sending data size = " << length << " to " << source);
  //buffer = new int[length];
  NEW_INT_ARRAY(buffer, length, *logger);
  buffer[0] = RT_NEIGH_DATA;
  buffer++; //skip RT_NEIGH_DATA
  prepare_neighbor_data_to_send(source, global_vid, buffer, length - 1);
  buffer--;

  DEBUG(*logger, "sending data");
  send_msg(buffer, length, source);
  DEBUG(*logger, "rank " << rank << " sent data size = " << length << " to " << source << " for global id  " << global_vid);
  //delete [] buffer;
  DELETE_INT_ARRAY(buffer, *logger);
}


void external_neighbor_handler::process_neighbor_request(int source, int *req, int req_len)
{
  DEBUG(*logger, "Got Neighbor request of size " << req_len << " from " << source);

  req++; req_len--; //skip RT_NEIGH_REQ
  int msg_count = 2;
  int total_messages = req_len / msg_count;
  if(total_messages * msg_count < req_len) total_messages++;

  for(int k = 0; k < total_messages; k++ ) {
    int *buffer = 0;
    std::vector<int> lengths;
    int total_length = 2;       //RT_NEIGH_DATA, num requests

    for(int i = k * msg_count; i < (k + 1) * msg_count && i< req_len; i++) {
      int global_vid = req[i];
      int degree = get_num_neighbors(source, global_vid);
      lengths.push_back(2 + degree * 4);
      total_length += lengths[i - k * msg_count];
    }

    buffer = new int[total_length];
    //NEW_INT_ARRAY(buffer, total_length, *logger);
    buffer[0] = RT_NEIGH_DATA;
    if ( (k + 1) * msg_count <= req_len)
      buffer[1] = msg_count;             //num_requests
    else
      buffer[1] = req_len - k * msg_count;           //num requests

    buffer += 2;       //skip RT_NEIGH_DATA , num requests
    DEBUG(*logger, "rank " << rank << " sending data size = " << total_length << " to " << source);

    stringstream ss;
    ss << "SENDING global vids (neighbor count)  to " << source << " : ";
    int offset = 0;
    for(int i = k * msg_count; i < (k + 1) * msg_count && i< req_len; i++) {
      int global_vid = req[i];
      prepare_neighbor_data_to_send(source, global_vid, (buffer + offset), lengths[i - k * msg_count]);
      ss << global_vid << "(" << ( (lengths[i - k * msg_count] - 2) / 4) << ") ";

      /*if (global_vid == 1434){
              cout<<"**** global vid 1434 ==> ";

              for(int l= offset; l < offset +lengths[i - k * msg_count] ; l++)
                      cout << buffer[l] <<" ";
              cout<<endl;
         }

         if (global_vid == 1368){
              cout<<"**** global vid 1368 ==> ";

              for(int l= offset; l < offset +lengths[i - k * msg_count] ; l++)
                      cout << buffer[l] <<" ";
              cout<<endl;
         }*/

      offset += lengths[i - k * msg_count];
    }

    //if (rank == 3)
    //INFO(*logger, ss.str());

    buffer -= 2;
    DEBUG(*logger, "sending data");
    DEBUG(*logger, "rank " << rank << " sent data size = " << total_length << " to " << source);
    send_msg(buffer, total_length, source);
    delete [] buffer;
    //DELETE_INT_ARRAY(buffer, *logger);
    //break;
  }

}

void external_neighbor_handler::process_request(int source, int size){

  int *recv_buf = new int[size];
  //int *recv_buf;
  //NEW_INT_ARRAY(recv_buf, size, *logger);
  int global_vid;
  MPI_Status stat;
  recv_msg(recv_buf, size, source, stat);
  DEBUG(*logger, "got message : " << get_request_type((REQUEST_TYPE)recv_buf[0])
        << ", from: " << source
        << "(=" << stat.MPI_SOURCE << ")");
  switch(recv_buf[0]) {
  case RT_NEIGH_REQ:
    /*global_vid = recv_buf[1];
       DEBUG(*logger, "rank " <<rank << " got neighbor request from: " << source << " for global id = " << global_vid);
       process_neighbor_request(source, global_vid);*/
    process_neighbor_request(source, recv_buf, size);
    break;
  case RT_NEIGH_DATA:
    TRACE4(*logger, "got RT_NEIGH_DATA");
    recv_buf++; //skip RT_NEIGH_DATA
    receive_neighbor_data(source, recv_buf, size - 1);
    recv_buf--;
    break;
  case RT_NEIGH_TOKEN:
    TRACE4(*logger, "got RT_NEIGH_TOKEN from " << source);
    if(rank == 0) {
      computation_finished = true;
      send_fin_to_all();
      has_token = false;   // to prevent sending fin the second time
    } else {
      has_token = true;
    }
    break;
  case RT_NEIGH_FIN:
    TRACE4(*logger, "got RT_NEIGH_FIN from " << source);
    computation_finished = true;
    break;
  default:
    CRITICAL_ERROR(*logger, "got an unknown request");
    delete[] recv_buf;
    exit(1);
  }

  delete[] recv_buf;
  //DELETE_INT_ARRAY(recv_buf, *logger);
}

void external_neighbor_handler::insert_pseudo_local(int local_id){
  pseudo_locals_on_extensions.insert(local_id);
}

void external_neighbor_handler::insert_pseudo_local_per_level(int level, int global_id){
  DEBUG(*logger, " global id " << global_id << " size " << pseudo_locals_on_extensions_per_level[level].size());
  omp_set_lock(&lock_nr);
  pseudo_locals_on_extensions_per_level[level].insert(global_id);
  omp_unset_lock(&lock_nr);
}

void external_neighbor_handler::handle_external_neighbors()
{

  send_neighbor_requests();
  external_neighbor_requests.clear();

  //omp_set_lock(&lock_nr);
  int src = -1;
  int size;
  do {
    src = check_message(size);
    if(src != -1) {
      process_request(src, size);
      //print_neighbor_request_list();
    }
  } while(!all_fetched());
  //reset_wait_for_external_neighbors();
  //omp_unset_lock(&lock_nr);
}

void external_neighbor_handler::handle_external_neighbors(int requester)
{
  DEBUG(*logger, "External neighbors to be requested = " << external_neighbor_requests.size());
  //return;
  if(numtasks > num_partitions)
    union_neighbor_requests_from_rank_group();

  computation_finished = false;
  if(rank == 0)
    has_token = true;
  else
    has_token = false;

  int src = -1;
  ranks_finished = 0;
  send_neighbor_requests();
  bool sent_fin = false;
  neighbors_fetched = 0;
  int size;
  do {
    src = check_message(size);
    if(src != -1) {
      process_request(src, size);
      //print_neighbor_request_list();
    }

    if(has_token && neighbors_fetched == external_neighbor_requests.size())
      forward_token();

  } while(computation_finished == false);

  DEBUG(*logger, "rank = " << rank << " FIN = " << ranks_finished);
  external_neighbor_requests.clear();
}

void external_neighbor_handler::add_neighbor_request(int part_id, int global_id){
  omp_set_lock(&lock_nr);
  external_neighbor_requests[global_id] = part_id;
  DEBUG(*logger, "Adding request for vertex " << global_id);
  omp_unset_lock(&lock_nr);
}

void external_neighbor_handler::add_neighbor_request(types::DFSCode dfscode, int part_id, int global_id){
  omp_set_lock(&lock_nr);
  external_neighbor_requests[global_id] = part_id;
  string code_str = dfscode.to_string();
  if(requests_for_dfscodes.count(code_str) == 0) {
    std::set<int> t;
    requests_for_dfscodes[code_str] = t;
  }
  requests_for_dfscodes[code_str].insert(global_id);

  DEBUG(*logger, "Adding request for vertex " << global_id);

  omp_unset_lock(&lock_nr);
}

void external_neighbor_handler::add_neighbor_request_no_lock(types::DFSCode dfscode, int part_id, int global_id){

  int thread_id = omp_get_thread_num();
  omp_set_lock(&lock_nr);
  external_neighbor_requests[global_id] = part_id;
  omp_unset_lock(&lock_nr);
  string code_str = dfscode.to_string();
  if(requests_for_dfscodes_threads[thread_id].count(code_str) == 0) {
    std::set<int> t;
    requests_for_dfscodes_threads[thread_id][code_str] = t;
  }
  requests_for_dfscodes_threads[thread_id][code_str].insert(global_id);

  DEBUG(*logger, "Adding request for vertex " << global_id);
}

void external_neighbor_handler::ext_handler_combine_threads_map_vectors(){
  for(int i = 0; i <numthreads; i++) {
    requests_for_dfscodes.insert(requests_for_dfscodes_threads[i].begin(), requests_for_dfscodes_threads[i].end());
    requests_for_dfscodes_threads[i].clear();
  }
}


void external_neighbor_handler::remove_neighbor_requests_for_non_frequent_dfscodes(std::vector<types::DFSCode> &dfscodes){

  std::map<int, int> filtered_requests;
  for(int i = 0; i < dfscodes.size(); i++) {
    string dfscode = dfscodes[i].to_string();
    if(requests_for_dfscodes.count(dfscode) > 0 ) {
      for(std::set<int>::iterator it = requests_for_dfscodes[dfscode].begin();
          it != requests_for_dfscodes[dfscode].end(); ++it) {
        filtered_requests[*it] = 0;                         // actual values are assigned in next step
      }
    }
  }

  for(std::map<int,int>::iterator it = filtered_requests.begin(); it != filtered_requests.end(); ++it) {
    it->second = external_neighbor_requests[it->first];
  }
  external_neighbor_requests.clear();
  external_neighbor_requests = filtered_requests;
}

//only thread 0 removes from the list, no locking needed
void external_neighbor_handler::remove_neighbor_request(int global_id){
  //omp_set_lock(&lock_nr);
  //std::set<int>::iterator it = external_neighbor_requests.find(global_id);
  //external_neighbor_requests.erase (it, external_neighbor_requests.end());
  //is_external_neighbor_fetched[global_id] = true;
  is_external_neighbor_fetched.erase(global_id);
  //omp_unset_lock(&lock_nr);
}

bool external_neighbor_handler::is_fetched(int global_id){
  omp_set_lock(&lock_nr);
  bool flag;
  if(is_external_neighbor_fetched.count(global_id) == 0)
    flag = false;
  flag = is_external_neighbor_fetched[global_id];
  omp_unset_lock(&lock_nr);
  return flag;
}

bool external_neighbor_handler::is_fetched(std::set<int> global_ids){

  omp_set_lock(&lock_nr);
  bool flag = true;
  for(std::set<int>::iterator it = global_ids.begin(); it != global_ids.end(); ++it) {
    if(is_external_neighbor_fetched.count(*it) == 0 || is_external_neighbor_fetched[*it] == false) {
      flag = false;
      break;
    }
  }
  omp_unset_lock(&lock_nr);
  return flag;
}

void external_neighbor_handler::print_neighbor_request_list(){
  omp_set_lock(&lock_nr);
  std::stringstream ss;
  ss << "rank " << rank << " waiting for neighbor requests : ";
  for(std::map<int,bool>::iterator it = is_external_neighbor_fetched.begin(); it != is_external_neighbor_fetched.end(); ++it) {
    if(it->second == false) {
      ss << " " << it->first;
    }
  }
  INFO(*logger, ss.str());
  omp_unset_lock(&lock_nr);
}

bool external_neighbor_handler::all_fetched(){
  omp_set_lock(&lock_nr);
  bool flag = true;
  for(std::map<int,bool>::iterator it = is_external_neighbor_fetched.begin(); it != is_external_neighbor_fetched.end(); ++it) {
    if(it->second == false) {
      flag = false;
      break;
    }
  }
  omp_unset_lock(&lock_nr);
  return flag;
}

void external_neighbor_handler::forward_token() {
  int buf[2];
  buf[0] = RT_NEIGH_TOKEN;
  buf[1] = 0;       //filler
  int dest = (rank + 1) % numtasks;
  send_msg(buf, 2, dest);
  has_token = false;
}

void external_neighbor_handler::send_fin_to_all() {
  int buf[2];
  buf[0] = RT_NEIGH_FIN;
  buf[1] = 0;       //filler
  for(int i = 1; i < numtasks; i++)
    send_msg(buf, 2, i);
}

void external_neighbor_handler::send_msg(int *buffer, int length, int dest_proc) {
  //omp_set_lock(&lock_nr);
  DEBUG(*logger, "XXX SENDING buffer of length " << length << " of type MSG_NEIGH_DATA to " << dest_proc);
  MPI_Isend(buffer, length, MPI_INT, dest_proc, MSG_NEIGH_DATA, MPI_COMM_WORLD, &request);
  //omp_unset_lock(&lock_nr);
}

void external_neighbor_handler::recv_msg(int *buffer, int length, int src_proc, MPI_Status &stat) {
  //omp_set_lock(&lock_nr);
  DEBUG(*logger, "XXX RECEIVING buffer of length " << length << " of type MSG_NEIGH_DATA from " << src_proc);
  MPI_Recv(buffer, length, MPI_INT, src_proc, MSG_NEIGH_DATA, MPI_COMM_WORLD, &stat);
  //omp_unset_lock(&lock_nr);
}

void external_neighbor_handler::all_gather_msg(int *sendbuf, int *recvbuf, int length, MPI_Comm rank_group) {
  DEBUG(*logger, dfscodes_sent << " XXX ALL_GATHER buffer of length " << length << " root " << root);
  MPI_Allgather(sendbuf, length, MPI_INT, recvbuf, length, MPI_INT, rank_group);
}

void external_neighbor_handler::all_gatherv_msg(int *sendbuf, int *recvbuf, int sendlen, int *recvlen, MPI_Comm rank_group) {
  DEBUG(*logger, dfscodes_sent << " XXX ALL_GATHERV buffer of length " << length << " root " << root);
  int *displs = new int[numtasks];
  displs[0] = 0;
  for(int i = 1; i< numtasks; i++)
    displs[i] = displs[i - 1] + recvlen[i - 1];
  MPI_Allgatherv(sendbuf, sendlen, MPI_INT, recvbuf, recvlen, displs, MPI_INT, rank_group);
  delete [] displs;
}

void external_neighbor_handler::alltoall_msg(int *sendbuf, int *recvbuf, int length, MPI_Comm rank_group) {
  DEBUG(*logger, dfscodes_sent << " XXX ALLTOALL buffer of length " << length << " root " << root);
  MPI_Alltoall(sendbuf, length, MPI_INT, recvbuf, length, MPI_INT, rank_group);
}


void external_neighbor_handler::alltoallv_msg(int *sendbuf, int *recvbuf, int *sendlen, int *recvlen, MPI_Comm rank_group) {
  DEBUG(*logger, dfscodes_sent << " XXX ALLTOALLV buffer of length " << length << " root " << root);
  int *sdispls = new int[numtasks];
  int *rdispls = new int[numtasks];
  sdispls[0] = 0;
  rdispls[0] = 0;
  for(int i = 1; i< numtasks; i++) {
    sdispls[i] = sdispls[i - 1] + sendlen[i - 1];
    rdispls[i] = rdispls[i - 1] + recvlen[i - 1];
  }

  MPI_Alltoallv(sendbuf, sendlen, sdispls, MPI_INT, recvbuf, recvlen, rdispls, MPI_INT, rank_group);
  delete [] rdispls;
  delete [] sdispls;
}


} //namespace alg
