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


#include <dynamic_hybrid_load_balancing_priv_params.hpp>
#include <unistd.h>

using namespace std;

namespace algs {

void dynamic_load_balancing::init_dynamic_load_balancing(unsigned int rank, unsigned int numtasks, MPI_Comm comm_object)
{
  this->rank = rank;
  this->numtasks = numtasks;
  //this->dfs_task_queue = dfs_task_queue;

  this->work_counter = 0;
  this->load_balance_interval = 5;
  this->do_process_load_balance = true;

  my_color = WHITE;
  has_token = false;
  token_color = WHITE;
  //next_work_request;
  computation_end = false;
  first_round = true;
  scheme = ARR;

  next_work_request = (rank + 1) % numtasks;
  if(rank == 0) has_token = true;

  //communication_object = comm_object.Dup();
  communication_object = comm_object;

  logger = Logger::get_logger("DYN_LOAD_BALANCE");
  TRACE(*logger, "dynamic load-balancing initialized");
}

void dynamic_load_balancing::reinitialize_dynamic_load_balancing()
{
  this->work_counter = 0;

  my_color = WHITE;
  has_token = false;
  token_color = WHITE;
  //next_work_request;
  computation_end = false;
  first_round = true;

  next_work_request = (rank + 1) % numtasks;
  if(rank == 0) has_token = true;

  TRACE(*logger, "dynamic load-balancing REinitialized");
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

void dynamic_load_balancing::set_process_load_balance_flag(bool flag){
  do_process_load_balance = flag;
}

bool dynamic_load_balancing::is_process_load_balance_on(){
  return do_process_load_balance;
}

void dynamic_load_balancing::send_work_request(){

  DEBUG(*logger, "trying to request work from: " << next_work_request << ", requested_work_from.size()=" << requested_work_from.empty());
  if(!requested_work_from.empty()) return;
  DEBUG(*logger, "is working = " << working() << " requesting work from: " << next_work_request);
  int buffer[2];
  buffer[0] = RT_WORK_REQUEST;
  send_other_msg(buffer, 2, next_work_request);
  requested_work_from.push_back(next_work_request);
}

/*
 * returns: -1 if no request is made, otherwise rank of the requesting processor
 */
int dynamic_load_balancing::check_request()
{
  MPI_Status stat;
  bool ret;
  TRACE5(*logger, "checking requests");
  if(computation_end == true) return -1;
  int flag;
  //MPI_Iprobe(MPI_ANY_SOURCE, MSG_REQUEST, MPI_COMM_WORLD, &flag, &stat);
  MPI_Iprobe(MPI_ANY_SOURCE, MSG_REQUEST, communication_object, &flag, &stat);

  //DEBUG(*logger, "prob value =" << ret );
  if(flag) {
    //TRACE(*logger, "got request from: " << stat.Get_source());
    /*if((MSG_TAG)stat.Get_tag() == MSG_DATA) {
       CRITICAL_ERROR(*logger, "got strange request with tag: MSG_DATA");
       //exit(1);
       }*/
    //return stat.Get_source();
    return stat.MPI_SOURCE;
  } else {
    //TRACE(*logger, "did not get any request, computation_end=" << computation_end);
    return -1;
  }
}


bool dynamic_load_balancing::receive_data(int source, int size, Thread_private_data &gprv)
{

  DEBUG(*logger, "has token=" << has_token);
  DEBUG(*logger, "receiving data of size=" << size << ", from = " << source);

  if(size == 0) {
    DEBUG(*logger, "did not get any work");

    //usleep(1000); //wait for microseconds;

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

  MPI_Status stat;
  if(requested_work_from[0] != source || requested_work_from.size() != 1) {
    CRITICAL_ERROR(*logger, "got request from: " << source
                   << ", while last request was from: " << requested_work_from[0]
                   << " STRANGE - EXITING");
    //exit(1);
  }

  int *itset_buffer = new int[size];

  recv_data_msg(itset_buffer, size, source, stat);
  int tmp_length = 0;

  process_received_data(itset_buffer, size, gprv);

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

  delete[] itset_buffer;
  return true;
}


bool dynamic_load_balancing::receive_data_parts(int source, int size, Thread_private_data &gprv)
{

  DEBUG(*logger, "has token=" << has_token);
  DEBUG(*logger, "receiving data of size=" << size << ", from = " << source);

  MPI_Status stat;
  if(requested_work_from[0] != source || requested_work_from.size() != 1) {
    CRITICAL_ERROR(*logger, "got request from: " << source
                   << ", while last request was from: " << requested_work_from[0]
                   << " STRANGE - EXITING");
    //exit(1);
  }

  int *itset_buffer = new int[size];

  int offset = 0;
  int buf[2];
  while(offset < size) {
    recv_other_msg(buf, 2, source, stat);
    int part_size = buf[1];
    recv_data_msg(itset_buffer + offset, part_size, source, stat);
    offset += part_size;
  }

  int tmp_length = 0;

  process_received_data(itset_buffer, size, gprv);

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

  delete[] itset_buffer;
  return true;
}



void dynamic_load_balancing::process_work_split_request(int source, Thread_private_data &gprv)
{
  // if this processor does not have work, send 0 (length of itemset) to the requesting processor
  if(working() == false) {
    int buffer[2];
    buffer[0] = RT_DATA;
    buffer[1] = 0;
    //communication_object.Send(&buffer, 2, MPI::INT, source, MSG_OTHER);
    send_other_msg(buffer, 2, source);
    DEBUG(*logger, "cannot split stack, refusing split stack request");
    return;
  }

  initiate_global_split_work(source, gprv);

  //Dijkstra's token termination is handled in complete_global_split_work
  //if(source < rank) my_color = BLACK;
}

void dynamic_load_balancing::process_request(int source, Thread_private_data &gprv){

  int recv_buf[2];
  MPI_Status stat;
  //communication_object.Recv(recv_buf, 2, MPI::INT, source, MSG_REQUEST, stat);
  recv_other_msg(recv_buf, 2, source, stat);
  DEBUG(*logger, "got request: " << get_request_type((REQUEST_TYPE)recv_buf[0])
        << ", from: " << source
        << "(=" << stat.Get_source() << ")");
  switch(recv_buf[0]) {
  case RT_WORK_REQUEST:
    DEBUG(*logger, "rank " << rank << " got work request from: " << source);
    process_work_split_request(source, gprv);
    break;
  case RT_TOKEN:
    DEBUG(*logger, "rank " << rank << " got token: " << get_token_color((TOKEN_COLOR)recv_buf[1]) << " from: " << source);
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
    DEBUG(*logger, "rank" << rank << " got RT_DATA from " << source);
    receive_data(source, recv_buf[1], gprv);
    return;
  case RT_DATA_PART:
    DEBUG(*logger, "rank" << rank << " got RT_DATA_PART from " << source);
    receive_data_parts(source, recv_buf[1], gprv);
    return;
  default:
    CRITICAL_ERROR(*logger, "got an unknown request");
    exit(1);
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

  DEBUG(*logger, "forwarding token from rank " << rank << " color " << get_token_color(token_color));

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

void dynamic_load_balancing::load_balance(Thread_private_data &gprv)
{

  //if(!do_process_load_balance)
  //return;

  if(numtasks == 1)
    return;

  int src = check_request();
  if(src != -1)
    process_request(src, gprv);

  //if(working() == true)
  //	DEBUG(*logger, " Processor " << rank <<" is working");

  if(working() == false) {
    if(has_token == true)
      process_token();

    if(do_process_load_balance == true && computation_end == false)
      send_work_request();
  }
}

void dynamic_load_balancing::load_balance(Thread_private_data &gprv, bool send_req)
{

  //if(!do_process_load_balance)
  //return;

  if(numtasks == 1)
    return;

  int src = check_request();
  if(src != -1)
    process_request(src, gprv);

  //if(working() == true)
  //	DEBUG(*logger, " Processor " << rank <<" is working");

  if(working() == false) {
    if(has_token == true)
      process_token();

    if(do_process_load_balance == true && computation_end == false && send_req)
      send_work_request();
  }
}


void dynamic_load_balancing::send_other_msg(int *buffer, int length, int dest_proc) {
  DEBUG(*logger, "XXX SENDING buffer of length " << length << " of type MSG_REQUEST to " << dest_proc
        << "; of type: " << get_request_type((REQUEST_TYPE) buffer[0]));
  MPI_Isend(buffer, length, MPI_INT, dest_proc, MSG_REQUEST, communication_object, &request);
}

void dynamic_load_balancing::send_data_msg(int *buffer, int length, int dest_proc) {
  DEBUG(*logger, "XXX SENDING buffer of length " << length << " of type MSG_DATA to " << dest_proc);
  MPI_Isend(buffer, length, MPI_INT, dest_proc, MSG_DATA, communication_object, &request);
}

void dynamic_load_balancing::recv_other_msg(int *buffer, int length, int src_proc, MPI_Status &stat) {
  DEBUG(*logger, "XXX RECEIVING buffer of length " << length << " of type MSG_REQUEST from " << src_proc);
  MPI_Recv(buffer, length, MPI_INT, src_proc, MSG_REQUEST, communication_object, &stat);
}

void dynamic_load_balancing::recv_data_msg(int *buffer, int length, int src_proc, MPI_Status &stat) {
  DEBUG(*logger, "XXX RECEIVING buffer of length " << length << " of type MSG_DATA from " << src_proc);
  MPI_Recv(buffer, length, MPI_INT, src_proc, MSG_DATA, communication_object, &stat);
}

bool dynamic_load_balancing::dijkstra_is_computation_finished() {
  if(rank == 0 && first_round == false && has_token && token_color == WHITE && working() == false)
    return true;
  return false;
}

} //namespace alg
