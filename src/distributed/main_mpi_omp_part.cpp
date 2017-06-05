/*
 * This source is part of the single graph mining algorithm.
 *
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

#include <iostream>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#include <logger.hpp>
#include <graph_output.hpp>
#include <graph_types.hpp>
#include <utils.hpp>
#include <dbio.hpp>
#include <memory_checker.hpp>
#include <graph_miner_mpi_omp_part.hpp>

using std::cout;
using std::endl;
using std::cerr;
using types::Graph;

Logger *value_logger = Logger::get_logger("VAL");
Logger *logger = Logger::get_logger("MAIN");

int main(int argc, char **argv)
{
  if(argc != 9 && argc != 7 && argc != 5) {
    std::cerr << "usage: \n" << argv[0] << " <filetype> <graph-filename> <support-int-abs> <num-partitions> [-p <partition-filename>] [-o <output-filename>]" << std::endl;
    return 1;
  }

  string filetype = argv[1];
  string filename = argv[2];
  int absolute_min_support = atoi(argv[3]);
  int num_partitions = atoi(argv[4]);
  string partition_filename, output_filename;

  if(argc >= 7) {
    if (strcmp(argv[5],"-p") == 0) {
      partition_filename = argv[6];
    } else if ( strcmp(argv[5],"-o") == 0) {
      output_filename = argv[6];
    } else {
      std::cerr << "usage: \n" << argv[0] << " <filetype> <graph-filename> <support-int-abs> <num-partitions> [-p <partition-filename>] [-o <output-filename>]" << std::endl;
      return 1;
    }
  }
  if(argc == 9) {
    if (strcmp(argv[7],"-o") == 0) {
      output_filename = argv[8];
    } else {
      std::cerr << "usage: \n" << argv[0] << " <filetype> <graph-filename> <support-int-abs> <num-partitions> [-p <partition-filename>] [-o <output-filename>]" << std::endl;
      return 1;
    }
  }


  if(absolute_min_support < 1) {
    cerr << "error: absolute_min_support < 1" << endl;
    return 3;
  }

  DEBUG(*logger, " Num Threads = " << omp_get_max_threads() );

  dbio::FILE_TYPE file_format;
  file_format = dbio::ftype2str(string(argv[1]));

  if(file_format == dbio::ADJP && partition_filename.size() == 0) {
    std::cerr << "error: specify the partition-filename with ADJP format" << std::endl;
    return 3;
  }else if (file_format == dbio::TXT){
	std::cerr << "txt format is not supported, use adjp or ladjp format." <<std::endl;
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
  //MPI_Get_processor_name(hostname, &len);
  //printf ("Number of tasks= %d My rank= %d Running on %s\n", numtasks, rank, hostname);

  if(num_partitions < 1 || numtasks % num_partitions != 0) {
    std::cerr << "error: numtasks must be divisible by the number of partitions" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  int ranks_per_partition = numtasks / num_partitions;

  Logger::set_parallel_info( rank, numtasks, true);

  Graph graph;

  if(file_format == dbio::LADJP) {
    std::stringstream ss;
    ss << filename << "/" << (rank % num_partitions);
    filename = ss.str();
  }

  if(file_format == dbio::ADJP) {
    dbio::read_graph_adj_par(filename, partition_filename, rank % num_partitions, graph);
  }
  else {
    dbio::read_graph(file_format, filename, graph);
  }

  graph_counting_output *count_out = new graph_counting_output();
  graph_file_output *file_out = 0;

  graph_output *current_out = count_out;

  if(rank == 0 && output_filename.size() > 0 && output_filename != "-") {
    //char rank_str[10];
    //sprintf(rank_str, "%d", rank);
    //file_out = new graph_file_output(output_filename + "p" + string(rank_str));
    //current_out = file_out;
    current_out = new graph_file_output(output_filename);
  }

  timeval start_time, stop_time;

  gettimeofday(&start_time, 0);

  LOG_VAL(*value_logger, absolute_min_support);
  LOG_VALUE(*value_logger, "output", current_out->to_string());

  GRAPH_MINER::graph_miner_mpi_omp_part fsm(1, num_partitions);      //threads per rank

  int epsilon = ceil(absolute_min_support * 0.0);
  fsm.set_graph(graph);
  fsm.set_min_support(absolute_min_support + epsilon);
  //char time_logfilename[1000];
  //sprintf(time_logfilename, "time-log-s%d-p%d-r%d",absolute_min_support, numtasks, rank);
  //fsm.set_time_logger(time_logfilename);

  char vmap_fprefix[1000];
  sprintf(vmap_fprefix, "vmaps-s%d-p%d-r%d", absolute_min_support, numtasks, rank);
  fsm.set_vmap_filename_prefix(std::string(vmap_fprefix));
  fsm.set_graph_output(current_out);

  fsm.run();

  gettimeofday(&stop_time, 0);

  DEBUG(*logger," clean-up...");

  MPI_Barrier( MPI_COMM_WORLD);

  int root = 0;
  double sendtime = utils::get_time_diff(start_time, stop_time), recvtime = 0.0;
  MPI_Reduce( &sendtime, &recvtime, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD );
  if (rank == root) {
    LOG_VALUE(*value_logger, "total_time", recvtime);
  }

  int sendbuf = fsm.get_frequent_patterns_count(), recvbuf = 0;
  MPI_Reduce( &sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD );
  if (rank == root) {
    LOG_VALUE(*value_logger, "total_frequent_graphs", recvbuf);
  }

  delete count_out;
  delete file_out;

  DEBUG(*logger, "maximal device memory usage (MB): " << memory_checker::get_max_memory_usage_mb());
  DEBUG(*logger, "clean-up...");
  //memory_checker::detect_memory_leaks();

  MPI_Finalize();

  return 0;
}
