# README #

* This contains implementations of parallel frequent graph mining algorithms for a single input graph. It includes sequential, parallel and distributed algorithms (with partitioned input graph) 
* The details of the algorithms be found in the disseration [Parallel Subgraph Mining on Hybrid Platforms: HPC Systems, Multi-cores and GPUs](http://www.cs.rpi.edu/~zaki/PaperDir/PhdTheses/talukder-thesis.pdf)
* Contact: Nilothpal Talukder, Email: <nilothpal.talukder@gmail.com>, Mohammed J. Zaki, Email: <zaki@cs.rpi.edu>

# Compile and Run instructions #

## Enter the directory and build ##

    cd DistGraph
    make

## Run the sequential miner ##

    ./src/sequential/graph_miner_seq -txt ../../testdata/pdb1.stxt 200

## Run the parallel miners ##

    mpirun -np 4 ./src/parallel/pargraph_mpi_dyn -txt ../../testdata/pdb_graph_single2.stxt 200
    mpirun -np 4 ./src/parallel/pargraph_mpi_split_all -txt ../../testdata/pdb_graph_single2.stxt 200
   
    #specify number of OpenMP threads
    export OMP_NUM_THREADS=4
    mpirun -np 4 ./src/parallel/pargraph_hybrid -txt ../../testdata/pdb1.stxt 200    

<code> pargraph_mpi_dyn</code> and <code>pargraph_mpi_split_all</code> performs MPI based parallelism and dynamic load balancing. When a work request is made from an idle process to an active process, the local task queue is split using two different strategies. The former one (<code> pargraph_mpi_dyn</code>) shares work from the bottom of the expanded pattern lattice, while the latter one (<code>pargraph_mpi_split_all</code>) shares work from all DFS levels of the lattice.
The hybrid parallel miner algorithm <code>pargraph_hybrid</code> uses two levels of parallelism (MPI+OpenMP) and uses hybrid dynamic load balancing. 
For description of the algorithms, please refer to the chapter 4 of the dissertation. 

## Run the distributed miner with graph file (adjacency list format) and partition file as input ##

    mpirun -np 4 ./src/distributed/distgraph -adjp ../../testdata/pdb1.metis 200 4 -p ../../testdata/pdb1.metis.part.4

## Run the distributed miner with local partitions as input ##

    mpirun -np 4 ./src/distributed/distgraph -ladjp ../../testdata/pdb1_part_4 200 4

<code>distgraph</code> takes partitioned input graph and mines in a distributed fashion using two levels of parallelism (MPI+OpenMP).
For detailed description, please refer to the chapter 3 of the dissertation. 

##input graph formats (txt, adjp and ladjp) ##
### txt format###
t # 0   
v 0 0   
v 1 1    
v 2 0    
v 3 1    
v 4 1    
v 5 0    
v 6 0    
v 7 1    
e 0 1 0    
e 0 2 0    
e 1 2 0    
e 2 3 0    
e 2 6 0    
e 3 4 0    
e 3 5 0    
e 4 5 0    
e 4 7 0    
e 5 6 0    
e 5 7 0    
e 6 7 0    

### adjp a.k.a. METIS format ###
The file format is as follows.
First line is : |V| |E| 011.   
The following |V| lines contain information and adjacency list of vertices 1 to |V|.   
Each line contains: vertex_label to_vertex_id1 edge_label1 to_vertex_id2 edge_label2 ... 
Note that vertex_ids and labels start from 1.

8 12 011   
1 2 1 3 1   
2 1 1 3 1   
1 1 1 2 1 4 1 7 1   
2 3 1 5 1 6 1   
2 8 1 4 1 6 1   
1 8 1 4 1 5 1 7 1   
1 8 1 3 1 6 1   
2 5 1 6 1 7 1   

### Partition file ### 
Our <code>distgraph</code> algorithm requires disjoint vertex partitioning. We can provide a partition file as input.
Each line indicate the partition id for the corresponding vertex_id.
In the file, the partition ids 0 to K-1 are assigned to vertex_ids 1 to |V|.
The output partition file from METIS also follows the same format.
The following is an example of 2 partitions of the above input graph.

0   
0   
0   
0   
1   
1   
1   
1   

### ladjp format ###
This format requires local partitions (with one hop external edge overlap) as file input. 
The files must be in the same directory, one file per partition. The filenames 0 to K-1 are used for K partitions. File format is as follows.
First line: |V|  (includes both internal and external vertices of the initial partition)
The following |V| lines contain information and adjacency list of vertices 0 to |V|-1 (local ids)   
Each line contains: partition_id global_vertex_id vertex_label to_local_vertex_id1 edge_label1 to_local_vertex_id2 edge_label2 ... 
Note: vertex_ids and labels start from 0!

Filename: dir/0  

7    
0 0 0 1 0 2 0    
0 1 1 0 0 2 0    
0 2 0 0 0 1 0 3 0 4 0    
0 3 1 2 0 5 0 6 0    
1 6 0 2 0    
1 4 1 3 0    
1 5 0 3 0    

Filename: dir/1

6    
1 4 1 3 0 4 0 1 0    
1 5 0 3 0 4 0 0 0 2 0   
1 6 0 3 0 5 0 1 0   
1 7 1 0 0 1 0 2 0   
0 3 1 0 0 1 0   
0 2 0 2 0   
