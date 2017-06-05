/*
 * This source is part of the single graph mining algorithm.
 *
 * Copyright 2013 Robert Kessl
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

#include <dynamic_load_balancing.hpp>
#include <unistd.h>

using namespace std;

namespace algs {

void dynamic_load_balancing::init_dynamic_load_balancing(unsigned int rank, unsigned int numtasks, MPI::Intracomm comm_object)
{
  this->rank = rank;
  this->numtasks = numtasks;

  this->work_counter = 0;
  this->load_balance_interval = 5;

  my_color = WHITE;
  has_token = false;
  token_color = WHITE;
  //next_work_request;
  computation_end = false;
  first_round = true;
  scheme = ARR;

  next_work_request = (rank + 1) % numtasks;
  if(rank == 0) has_token = true;

  communication_object = comm_object.Dup();

  logger = Logger::get_logger("DYN_LOAD_BALANCE");
  TRACE(*logger, "dynamic load-balancing initialized");
}


void dynamic_load_balancing::set_load_balance_interval(int i){
  load_balance_interval = i;
}

bool dynamic_load_balancing::set_load_balancing_scheme(DYN_SCHEME scheme){
  if( scheme != ARR or scheme != RP)
    return false;
  this->scheme = scheme;

  if(scheme == RP) {
    next_work_request = random() % numtasks;
    while(next_work_request == rank)
      next_work_request = random() % numtasks;
  }
}

void dynamic_load_balancing::send_work_request(){

  DEBUG(*logger, "trying to request work from: " << next_work_request << ", requested_work_from.size()=" << requested_work_from.empty());
  if(!requested_work_from.empty()) return;
  TRACE(*logger, "is working = " << working() << " requesting work from: " << next_work_request);
  int buffer[2];
  buffer[0] = RT_WORK_REQUEST;
  //MPI::Status stat;
  //communication_object.Send(buffer, 2, MPI::INT, next_work_request, MSG_OTHER);
  send_other_msg(buffer, 2, next_work_request);
  requested_work_from.push_back(next_work_request);
}

/*
 * returns: -1 if no request is made, otherwise rank of the requesting processor
 */
int dynamic_load_balancing::check_request(bool blocking)
{
  MPI::Status stat;
  bool ret;
  TRACE5(*logger, "checking requests");
  if(computation_end == true) return -1;
  if(blocking == false)
    //ret = communication_object.Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, stat);
    ret = communication_object.Iprobe(MPI_ANY_SOURCE, MSG_REQUEST, stat);
  else {
    communication_object.Probe(MPI_ANY_SOURCE, MSG_REQUEST, stat);
    DEBUG(*logger, "got request from: " << stat.Get_source());
    /*if((MSG_TAG)stat.Get_tag() == MSG_DATA) {
       CRITICAL_ERROR(*logger, "got strange request with tag: MSG_DATA");
       //exit(1);
       }*/
    return stat.Get_source();
  }

  DEBUG(*logger, "prob value =" << ret );
  if(ret == true) {
    TRACE(*logger, "got request from: " << stat.Get_source());
    if((MSG_TAG)stat.Get_tag() == MSG_DATA) {
      CRITICAL_ERROR(*logger, "got strange request with tag: MSG_DATA");
      //exit(1);
    }
    return stat.Get_source();
  } else {
    //TRACE(*logger, "did not get any request, computation_end=" << computation_end);
    return -1;
  }
}


bool dynamic_load_balancing::receive_data(int source, int size)
{

  DEBUG(*logger, "has token=" << has_token);
  TRACE(*logger, "receiving data of size=" << size << ", from = " << source);

  if(size == 0) {
    TRACE(*logger, "did not get any work");

    usleep(1000); //wait for microseconds;

    if(scheme == RP) {
      next_work_request = random() % numtasks;
      while(next_work_request == rank)
        next_work_request = random() % numtasks;
    }else{ //ARR
      next_work_request = (next_work_request + 1) % numtasks;
      if(next_work_request == rank)  //make sure that the request is not made to self
        next_work_request = (next_work_request + 1) % numtasks;
    }

    requested_work_from.erase(requested_work_from.begin());

    return false;
  }

  //int buffer_size[2];
  MPI::Status stat;

  //if()
  if(requested_work_from[0] != source || requested_work_from.size() != 1) {
    CRITICAL_ERROR(*logger, "got request from: " << source
                   << ", while last request was from: " << requested_work_from[0]
                   << " STRANGE - EXITING");
    //exit(1);
  }

  int *itset_buffer = new int[size];

  //communication_object.Recv(itset_buffer, size, MPI::INT, source, MSG_DATA, stat);
  recv_data_msg(itset_buffer, size, source, stat);
  int tmp_length = 0;
  //communication_object.GetCount(&stat, MPI::INT, &tmp_length);
  //tmp_length = stat.Get_count(MPI::INT);
  //TRACE5(*logger, "MPI status says: " << tmp_length << " integers, request says: " << size << " integers");

  //deserialize_itemset(itset_buffer, FI_tree::items_to_compute);
  //TRACE(*logger, "got " << FI_tree::items_to_compute.size() << " items from #" << next_work_request);
  process_received_data(itset_buffer, size);
  if(scheme == RP) {
    next_work_request = random() % numtasks;
    while(next_work_request == rank)
      next_work_request = random() % numtasks;
  }else{ //ARR
    next_work_request = (next_work_request + 1) % numtasks;
    if(next_work_request == rank)  //make sure that the request is not made to self
      next_work_request = (next_work_request + 1) % numtasks;
  }

  requested_work_from.erase(requested_work_from.begin());

  //deleted in process_received_data
  //delete[] itset_buffer;

  return true;
}

void dynamic_load_balancing::process_work_split_request(int source)
{
  // if this processor does not have work, send 0 (length of itemset) to the requesting processor
  if(working() == false || can_split_work() == false) {
    int buffer[2];
    buffer[0] = RT_DATA;
    buffer[1] = 0;
    //communication_object.Send(&buffer, 2, MPI::INT, source, MSG_OTHER);
    send_other_msg(buffer, 2, source);
    DEBUG(*logger, "cannot split stack, refusing split stack request");
    return;
  }

  int *buffer = 0;
  int length = 0;
  split_work(buffer, length);

  int buffer_size[2];
  buffer_size[0] = RT_DATA;
  buffer_size[1] = length; // put there length of the split stack split.size()+1;

  TRACE(*logger, "sending data size = " << length << "to " << source);
  //communication_object.Send(&buffer_size, 2, MPI::INT, source, MSG_OTHER);
  send_other_msg(buffer_size, 2, source);

  DEBUG(*logger, "sending data");
  //communication_object.Send(buffer, split.size()+1, MPI::INT, source, MSG_DATA);
  //communication_object.Send(buffer, length, MPI::INT, source, MSG_DATA);
  send_data_msg(buffer, length, source);

  //handle Dijkstra's token termination
  if(source < rank) my_color = BLACK;
  delete [] buffer;

}

void dynamic_load_balancing::process_request(int source){

  int recv_buf[2];
  MPI::Status stat;
  //communication_object.Recv(recv_buf, 2, MPI::INT, source, MSG_REQUEST, stat);
  recv_other_msg(recv_buf, 2, source, stat);
  DEBUG(*logger, "got request: " << get_request_type((REQUEST_TYPE)recv_buf[0])
        << ", from: " << source
        << "(=" << stat.Get_source() << ")");
  switch(recv_buf[0]) {
  case RT_WORK_REQUEST:
    DEBUG(*logger, "got work request from: " << source);
    process_work_split_request(source);
    break;
  case RT_TOKEN:
    TRACE(*logger, "got token: " << get_token_color((TOKEN_COLOR)recv_buf[1]) << " from: " << source);
    token_color = (TOKEN_COLOR)recv_buf[1];
    has_token = true;
    break;
  case RT_COMPUTATION_END:
    TRACE(*logger, "got COMPUTATION_END");
    if(working())
      TRACE(*logger, "got COMPUTATION_END while working");
    computation_end = true;
    return;
  case RT_DATA:
    TRACE4(*logger, "got RT_DATA");
    receive_data(source, recv_buf[1]);
    return;
  default:
    CRITICAL_ERROR(*logger, "got an unknown request");
    //exit(1);
  }
}

/** if the processor has token follow the Dijkstra's computation end algorithm and forward a token with
 * the appropriate color
 */

void dynamic_load_balancing::process_token()
{

  if(dijkstra_is_computation_finished()) {       //only for rank 0
    //broadcast termination to everyone
    initiate_termination();
    return;
  }

  if(first_round == true)
    first_round = false;

  //for processor 0 always forward white token
  if(rank == 0)
    token_color = WHITE;

  TRACE(*logger, "forwarding token from rank " << rank << " color " << get_token_color(token_color));

  if(my_color == WHITE)
    forward_token(token_color);
  else
    forward_token(BLACK);
  has_token = false;
  my_color = WHITE;

}

void dynamic_load_balancing::forward_token(TOKEN_COLOR tcol) {
  int buffer[2];
  buffer[0] = RT_TOKEN;
  switch(tcol) {
  case BLACK:
    buffer[1] = BLACK;
    break;
  case WHITE:
    buffer[1] = WHITE;
    break;
  } // switch
  send_other_msg(buffer, 2, (rank + 1) % numtasks);
} // forward_token


void dynamic_load_balancing::initiate_termination(){

  TRACE(*logger, "computation end, sending COMPUTATION_END to all processors");
  //replace with broadcast
  for(int i = 1; i < numtasks; i++) {
    int buffer[2];
    buffer[0] = RT_COMPUTATION_END;
    //communication_object.Send(&buffer, 2, MPI::INT, i, MSG_OTHER);
    send_other_msg(buffer, 2, i);
  }
  computation_end = true;
}

void dynamic_load_balancing::load_balance()
{
  work_counter++;
  if(work_counter % load_balance_interval != 0)
    return;

  int src = check_request(false);
  if(src != -1)
    process_request(src);

  //if(working() == true)
  //	DEBUG(*logger, " Processor " << rank <<" is working");

  if(working() == false) {
    if(has_token == true)
      process_token();

    if(computation_end == false)
      send_work_request();
  }
}


void dynamic_load_balancing::send_other_msg(int *buffer, int length, int dest_proc) {
  TRACE5(*logger, "XXX SENDING buffer of length " << length << " of type MSG_REQUEST to " << dest_proc
         << "; of type: " << get_request_type((REQUEST_TYPE) buffer[0]));
  communication_object.Isend(buffer, length, MPI::INT, dest_proc, MSG_REQUEST);
}

void dynamic_load_balancing::send_data_msg(int *buffer, int length, int dest_proc) {
  TRACE5(*logger, "XXX SENDING buffer of length " << length << " of type MSG_DATA to " << dest_proc);
  communication_object.Isend(buffer, length, MPI::INT, dest_proc, MSG_DATA);
}

void dynamic_load_balancing::recv_other_msg(int *buffer, int length, int src_proc, MPI::Status &stat) {
  TRACE5(*logger, "XXX RECEIVING buffer of length " << length << " of type MSG_REQUEST from " << src_proc);
  communication_object.Recv(buffer, length, MPI::INT, src_proc, MSG_REQUEST, stat);
}

void dynamic_load_balancing::recv_data_msg(int *buffer, int length, int src_proc, MPI::Status &stat) {
  TRACE5(*logger, "XXX RECEIVING buffer of length " << length << " of type MSG_DATA from " << src_proc);
  communication_object.Recv(buffer, length, MPI::INT, src_proc, MSG_DATA, stat);
}

bool dynamic_load_balancing::dijkstra_is_computation_finished() {
  if(rank == 0 && first_round == false && has_token && token_color == WHITE && working() == false)
    return true;
  return false;
}


} //namespace alg
