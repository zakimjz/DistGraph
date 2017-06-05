export OMP_NUM_THREADS=4
srun --nodes=4 --ntasks-per-node=1 --cpus-per-task=4 ../src/distributed/distgraph -ladjp ../testdata/pdb1_part_4 200 4 -o output
