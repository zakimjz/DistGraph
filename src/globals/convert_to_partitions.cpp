#include <string>
#include <iomanip>
#include <iostream>

#include <graph_output.hpp>
#include <graph_types.hpp>
#include <logger.hpp>
#include <utils.hpp>
#include <dbio.hpp>
#include <mpi.h>
#include <math.h>

using std::string;
using namespace types;
using namespace dbio;
using std::cerr;
using std::fixed;

Logger *value_logger = Logger::get_logger("VAL");
Logger *logger = Logger::get_logger("MAIN");

int main(int argc, char **argv)
{
  if(argc != 5) {
    std::cerr << "usage: \n" << argv[0] << " <adj-filename-in>  <partition-filename-in> <num-partitions> <dirname-out>" << std::endl;
    return 1;
  }

  string graph_partition_in = argv[1];
  string partition_filename_in = argv[2];
  int num_partitions = atoi(argv[3]);
  string dirname_out = argv[4];

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

  for(int i = 0; i < num_partitions; i++){
	  if( i % numtasks == rank ){
		  Graph graph;
		  std::cout << "reading graph rank = " << rank << " partition id = " << i << endl;
		  dbio::read_graph_adj_par(graph_partition_in, partition_filename_in, i, graph);
		  std::cout << "writing graph rank = " << rank << " partition id = " << i << endl;
		  dbio::write_local_adjp(dirname_out, i, graph);
	  }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

  return 0;
}
