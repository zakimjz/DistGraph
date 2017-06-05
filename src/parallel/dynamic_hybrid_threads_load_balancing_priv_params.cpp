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


#include <dynamic_hybrid_threads_load_balancing_priv_params.hpp>
#include <unistd.h>

using namespace std;

namespace algs {

void dynamic_threads_load_balancing::init_dynamic_threads_load_balancing( int rank, int num_threads)
{
  this->rank = rank;
  this->num_threads = num_threads;
  this->scheme = RP;

  next_work_request.clear();
  message_queue.clear();
  requested_work_from.clear();

  for(int i = 0; i < num_threads; i++) {
    int p =  (i + 1) % num_threads;
    next_work_request.push_back(p);

    std::deque<int> dq;
    message_queue.push_back(dq);
    std::vector<int> vc;
    requested_work_from.push_back(vc);
  }
  all_idle_work_request_counter = 0;
  logger = Logger::get_logger("DYN_THREADS_LOAD_BALANCE");
  TRACE(*logger, "dynamic threads load-balancing initialized");
}

void dynamic_threads_load_balancing::set_num_threads(int num_threads){
  this->num_threads = num_threads;
}

void dynamic_threads_load_balancing::reinitialize_dynamic_threads_load_balancing()
{
  next_work_request.clear();
  for(int i = 0; i < num_threads; i++) {
    int p =  (i + 1) % num_threads;
    next_work_request.push_back(p);
  }
  all_idle_work_request_counter = 0;
  TRACE(*logger, "dynamic threads load-balancing REinitialized");
}

bool dynamic_threads_load_balancing::set_load_balancing_scheme(DYN_SCHEME scheme){
  if( scheme != ARR or scheme != RP)
    return false;
  this->scheme = scheme;

  int thread_id = omp_get_thread_num();
  if(scheme == RP) {
    next_work_request[thread_id] = random() % num_threads;
    while(next_work_request[thread_id] == thread_id)
      next_work_request[thread_id] = random() % num_threads;
  }
}

void dynamic_threads_load_balancing::send_work_request(Thread_private_data &gprv){

  int thread_id = gprv.thread_id; //omp_get_thread_num();
  if(!requested_work_from[thread_id].empty()) return;

  //omp_set_lock(&lock);
  DEBUG(*logger, "rank " << rank << " thread " << thread_id << " requesting work from thread " << next_work_request[thread_id]
        << ", requested_work_from.size() = " << requested_work_from[thread_id].size());
  //omp_unset_lock(&lock);

  int buffer[2];
  buffer[0] = RT_WORK_REQUEST;
  buffer[1] = 0;       // filler
  send_msg(buffer, 2, thread_id, next_work_request[thread_id]);
  requested_work_from[thread_id].push_back(next_work_request[thread_id]);
}

void dynamic_threads_load_balancing::send_work_request(Thread_private_data &gprv, int dest){

  int thread_id = gprv.thread_id; //omp_get_thread_num();
  if(!requested_work_from[thread_id].empty()) return;
  next_work_request[thread_id] =  dest;

  omp_set_lock(&lock);
  DEBUG(*logger, "thread " << thread_id << " IDLE: just send work request to 0");
  DEBUG(*logger, "rank " << rank << " IDLE; thread " << thread_id << " requesting work from thread " << next_work_request[thread_id]
        << ", requested_work_from.size() = " << requested_work_from[thread_id].size());
  omp_unset_lock(&lock);

  int buffer[2];
  buffer[0] = RT_WORK_REQUEST;
  buffer[1] = 0; // filler
  send_msg(buffer, 2, thread_id, dest);
  requested_work_from[thread_id].push_back(dest);
}

/*
 * returns: -1 if no request is made, otherwise rank of the requesting processor
 */
int dynamic_threads_load_balancing::check_request(Thread_private_data &gprv)
{
  int thread_id = gprv.thread_id; //omp_get_thread_num();
  TRACE5(*logger, "checking requests");
  if(message_queue[thread_id].size() > 0 ) {
    omp_set_lock(&lock);
    int source = message_queue[thread_id].front();
    DEBUG(*logger, "thread " << thread_id << " threads dyn load balancing, got request from: " << source);
    omp_unset_lock(&lock);
    return source;
  }
  return -1;
}


bool dynamic_threads_load_balancing::receive_data(int source, int size, Thread_private_data &gprv)
{
  int thread_id = gprv.thread_id; //omp_get_thread_num();

  //omp_set_lock(&lock);
  DEBUG(*logger, "rank " << rank << " thread " << thread_id << " receiving data of size=" << size << ", from thread = " << source);
  //omp_unset_lock(&lock);

  if(size == 0) {
    DEBUG(*logger, "thread" << thread_id << " did not get any work");
    //usleep(10000); //wait for microseconds;
    if(scheme == RP) {
      next_work_request[thread_id] = random() % num_threads;
      while(next_work_request[thread_id] == thread_id)
        next_work_request[thread_id] = random() % num_threads;
    }else{ //ARR
      next_work_request[thread_id] = (next_work_request[thread_id] + 1) % num_threads;
      if(next_work_request[thread_id] == thread_id) //make sure that the request is not made to self
        next_work_request[thread_id] = (next_work_request[thread_id] + 1) % num_threads;
    }

    requested_work_from[thread_id].erase(requested_work_from[thread_id].begin());
    return false;
  }

  /*if(requested_work_from[thread_id].size() == 0 && source != 0 ){
            CRITICAL_ERROR(*logger, "thread" << thread_id <<" got response from: " << source
                           << ", while last request was from: " << requested_work_from[thread_id][0]
                                           << " size: " << requested_work_from[thread_id].size()
                           << " STRANGE - EXITING");
            exit(1);
     }*/
  if(requested_work_from[thread_id].size() != 1 || requested_work_from[thread_id][0] != source ) {
    CRITICAL_ERROR(*logger, "thread" << thread_id << " got response from: " << source
                   << ", while last request was from: " << requested_work_from[thread_id][0]
                   << " size: " << requested_work_from[thread_id].size()
                   << " STRANGE - EXITING");
    exit(1);
  }

  ////old comment: nothing else to do, data is already pushed in the queue by the donor thread
  // process the data put in the shared queue
  thread_process_received_data(gprv);

  if(scheme == RP) {
    next_work_request[thread_id] = random() % num_threads;
    while(next_work_request[thread_id] == thread_id)
      next_work_request[thread_id] = random() % num_threads;
  }else{ //ARR
    next_work_request[thread_id] = (next_work_request[thread_id] + 1) % num_threads;
    if(next_work_request[thread_id] == thread_id) //make sure that the request is not made to self
      next_work_request[thread_id] = (next_work_request[thread_id] + 1) % num_threads;
  }

  requested_work_from[thread_id].erase(requested_work_from[thread_id].begin());

  thread_start_working();

  return true;
}

void dynamic_threads_load_balancing::process_work_split_request(int source, Thread_private_data &gprv)
{
  int thread_id = gprv.thread_id; //omp_get_thread_num();

  //if threadid = 0 and all threads are idle, no need to send response now
  // if global work is received the response will be sent eventually
  if(thread_id == 0 && all_threads_idle()) {
    all_idle_work_request_counter++;
    DEBUG(*logger, " thread " << source << " requested work in all IDLE state: ignored, request counter " << all_idle_work_request_counter);
    return;
  }

  if(thread_working(gprv) == false || can_thread_split_work(gprv) == false) {
    int buffer[2];
    buffer[0] = RT_WORK_RESPONSE;
    buffer[1] = 0;
    send_msg(buffer, 2, thread_id, source);
    DEBUG(*logger, "thread " << thread_id << " cannot split stack, refusing split stack request from thread " << source);

    /*
       next_work_request[source] = (next_work_request[source]+1) % num_threads;
       if(next_work_request[source] == source) //make sure that the request is not made to self
            next_work_request[source] = (next_work_request[source]+1)% num_threads;
       requested_work_from[source].erase(requested_work_from[source].begin());
     */

    return;
  }

  int length;
  thread_split_work(source, length, gprv);

  int buffer_size[2];
  buffer_size[0] = RT_WORK_RESPONSE;
  buffer_size[1] = length; // put there length of the split stack split.size()+1;

  DEBUG(*logger, "sending data size = " << length << "to " << source);
  send_msg(buffer_size, 2, thread_id, source);

  /*
     next_work_request[source] = (next_work_request[source]+1) % num_threads;
     if(next_work_request[source] == source) //make sure that the request is not made to self
            next_work_request[source] = (next_work_request[source]+1)% num_threads;
     requested_work_from[source].erase(requested_work_from[source].begin());
   */
}

void dynamic_threads_load_balancing::process_request(int source, Thread_private_data &gprv){

  int thread_id = gprv.thread_id; //omp_get_thread_num();

  int recv_buf[2];
  recv_msg(recv_buf, 2, thread_id, source);
  DEBUG(*logger, "thread " << thread_id << " got request: " << get_request_type((REQUEST_TYPE)recv_buf[0])
        << ", from: " << source );
  switch(recv_buf[0]) {
  case RT_WORK_REQUEST:
    DEBUG(*logger, "thread " << thread_id << " got work request from: " << source);
    process_work_split_request(source, gprv);
    break;
  case RT_WORK_RESPONSE:
    DEBUG(*logger, "thread " << thread_id << "got RT_WORK_RESPONSE");
    receive_data(source, recv_buf[1], gprv);
    return;
  case GLOBAL_WORK_REQUEST:
    DEBUG(*logger, "thread " << thread_id << " got global work request from: " << source);
    process_global_work_split_request(recv_buf[1], gprv);
    break;
  case GLOBAL_WORK_RESPONSE:
    DEBUG(*logger, "thread " << thread_id << "got GLOBAL_WORK_RESPONSE");
    process_global_work_split_response(recv_buf[1], gprv);
    break;
  default:
    CRITICAL_ERROR(*logger, "got an unknown request");
    exit(1);
  }

}

void dynamic_threads_load_balancing::threads_load_balance(Thread_private_data &gprv)
{

  int thread_id = gprv.thread_id; //omp_get_thread_num();
  int src = check_request(gprv);

  if(src != -1)
    process_request(src, gprv);

  if(thread_working(gprv) == false) {
    if(all_threads_idle()) {
      send_work_request(gprv, 0);
    }
    else
      send_work_request(gprv);            //global load balance
  }
}


void dynamic_threads_load_balancing::send_msg(int *buffer, int length, int src_thr, int dest_thr){

  omp_set_lock(&lock);
  DEBUG(*logger, "thread " << src_thr << " SENDING buffer of length " << length << " of type MSG_REQUEST to thread " << dest_thr
        << "; of type: " << get_request_type((REQUEST_TYPE) buffer[0]));
  message_queue[dest_thr].push_back(src_thr);
  for(int i = 0; i <length; i++)
    message_queue[dest_thr].push_back(buffer[i]);

  omp_unset_lock(&lock);
}

void dynamic_threads_load_balancing::recv_msg(int *buffer, int length, int thr, int originating_thr) {

  omp_set_lock(&lock);

  int source = message_queue[thr].front();
  if(originating_thr != source) {
    CRITICAL_ERROR(*logger, "thread " << thr << ": message sources did not match " << source << " and " << originating_thr);
    exit(0);
  }
  message_queue[thr].pop_front(); //take off source
  for(int i = 0; i < length; i++) {
    buffer[i] = message_queue[thr].front();
    message_queue[thr].pop_front();
  }

  DEBUG(*logger, "thread " << thr << " RECEIVING buffer of length " << length << " of type MSG_REQUEST from " << originating_thr
        << "; of type: " << get_request_type((REQUEST_TYPE) buffer[0]));

  omp_unset_lock(&lock);
}

/************************************************
*
* Hybrid Load Balancing Functions starts here
*
************************************************/

void dynamic_threads_load_balancing::send_global_split_requests(int requester_rank_id, Thread_private_data &gprv){
  //send all thread global split request
  int buf[2];
  buf[0] = GLOBAL_WORK_REQUEST;
  buf[1] = requester_rank_id;         // put there length of the split stack split.size()+1;

  for(int i = 0; i< num_threads; i++) {
    if(i != gprv.thread_id)
      send_msg(buf, 2, gprv.thread_id, i);
  }
  DEBUG(*logger, "thread " << gprv.thread_id << " sent global split requests");
}

void dynamic_threads_load_balancing::process_global_work_split_request(int requester_rank_id, Thread_private_data &gprv){
  thread_split_global_work(requester_rank_id, gprv);
  //send global work response to thread 0
  int buf[2];
  buf[0] = GLOBAL_WORK_RESPONSE;
  buf[1] = requester_rank_id;
  send_msg(buf, 2, gprv.thread_id, 0);
  DEBUG(*logger, "thread " << gprv.thread_id << " sent global split reply");
}

void dynamic_threads_load_balancing::process_global_work_split_response(int requester_rank_id, Thread_private_data &gprv){
  if(global_work_response_counter.count(requester_rank_id) == 0)
    global_work_response_counter[requester_rank_id] = 0;
  global_work_response_counter[requester_rank_id]++;
  if(global_work_response_counter[requester_rank_id] == num_threads - 1) {      //exclude itself (thread 0)
    complete_global_work_split_request(requester_rank_id, gprv);
    global_work_response_counter[requester_rank_id] = 0;
  }
}

void dynamic_threads_load_balancing::send_work_response_to_all_threads(Thread_private_data &gprv, std::vector<int> thread_has_work){
  //send all threads work response
  int buf[2];
  buf[0] = RT_WORK_RESPONSE;

  for(int i = 0; i< num_threads; i++) {
    if(i != gprv.thread_id) {            //0
      buf[1] = thread_has_work[i];             // length, anything but 0 is fine
      send_msg(buf, 2, gprv.thread_id, i);
    }
  }
  DEBUG(*logger, "thread " << gprv.thread_id << " sent global work response to all threads");
}

bool dynamic_threads_load_balancing::all_idle_work_requests_collected(){
  return (all_idle_work_request_counter == num_threads);
}

void dynamic_threads_load_balancing::reset_idle_work_requests_counter(){
  all_idle_work_request_counter = 0;
}

} //namespace alg
