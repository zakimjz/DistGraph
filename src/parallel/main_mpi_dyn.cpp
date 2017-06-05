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


#include <iostream>
#include <logger.hpp>
#include <graph_output.hpp>
#include <graph_types.hpp>
#include <utils.hpp>
#include <dbio.hpp>

#include <graph_miner_mpi_dyn.hpp>

#include <sys/time.h>
#include <mpi.h>

using std::cout;
using std::endl;
using std::cerr;
using types::Graph;

Logger *value_logger = Logger::get_logger("VAL");
Logger *logger = Logger::get_logger("MAIN");

int main(int argc, char **argv)
{
  if(argc != 5 && argc != 4) {
    std::cerr << "usage: \n" << argv[0] << " <filetype> <filename> <support-int-abs> [<output-filename>]" << std::endl;
    return 1;
  }

  string filetype = argv[1];
  string filename = argv[2];
  int absolute_min_support = atoi(argv[3]);
  string output_filename;
  if(argc == 5) output_filename = argv[4];


  if(absolute_min_support < 1) {
    cerr << "error: absolute_min_support < 1" << endl;
    return 3;
  }

  /* MPI code begins */

  int numtasks, rank, len, rc;
  char hostname[MPI_MAX_PROCESSOR_NAME];

  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    cerr << "Error starting MPI program. Terminating." << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


  Logger::set_parallel_info( rank, numtasks, true);

  Graph graph;

  dbio::FILE_TYPE file_format;
  file_format = dbio::ftype2str(string(argv[1]));
  dbio::read_graph(file_format, filename, graph);
  graph_counting_output *count_out = new graph_counting_output();
  graph_file_output *file_out = 0;

  graph_output *current_out = count_out;

  if(output_filename.size() > 0 && output_filename != "-") {
	//file_out = new graph_file_output(output_filename);
	char rank_str[10];
	sprintf(rank_str, "%d", rank);
	file_out = new graph_file_output(output_filename + "p" + string(rank_str));
    current_out = file_out;
  }
  timeval start_time, stop_time;

  gettimeofday(&start_time, 0);

  LOG_VAL(*value_logger, absolute_min_support);
  LOG_VALUE(*value_logger, "output", current_out->to_string());

  //if(rank == 0){
  GRAPH_MINER::graph_miner_mpi_dyn fsm;

  fsm.set_graph(graph);
  fsm.set_min_support(absolute_min_support);
  //fsm.create_neighbor_map();
  fsm.set_task_split_threshold(2);
  fsm.set_load_balance_interval(5);
  fsm.set_graph_output(current_out);

  fsm.run();

  gettimeofday(&stop_time, 0);

  LOG_VALUE(*value_logger, "total_time_for_rank_" << rank, std::fixed << utils::get_time_diff(start_time, stop_time));
  LOG_VALUE(*value_logger, "total_frequent_graphs_for_rank_" << rank, current_out->get_count());
  INFO(*logger," clean-up...");


  MPI_Barrier( MPI_COMM_WORLD);

  int root;
  double sendtime = utils::get_time_diff(start_time, stop_time), recvtime = 0.0;
  MPI_Reduce( &sendtime, &recvtime, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD );
  if (rank == root) {
    LOG_VALUE(*value_logger, "total_time", recvtime);
  }

  int sendbuf = current_out->get_count(), recvbuf = 0;
  MPI_Reduce( &sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD );
  if (rank == root) {
    LOG_VALUE(*value_logger, "total_frequent_graphs", recvbuf);
  }

  delete count_out;
  delete file_out;

  MPI_Finalize();

  return 0;
}

