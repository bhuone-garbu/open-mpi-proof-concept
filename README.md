# open-mpi-proof-concept

This is proof of concept I wrote when I was at college/uni using OpenMPI to parallelize some computation on a 2-d array.
https://www.open-mpi.org/

Basically, this is an an array relaxation that computes an average values based on surrounding neighbour values until a set  precision is reached.
The array is initially setup at the root processor (with rank = 0) and the rows of the array is given to the set number of processors. Array rows are scatter in such a way that a processor can work in isolation without needing to talk to
other processors during its average value calculation process. Once every processors finish working on their arrays, processors will send and request new rows to and from other processor that might affect their next calculation. This process
is repeated until the precision has been met. Finally, the new rows are gathered
back the main array in the root processor.

When running the program with mpirun I have also added further two arguments to specify the dimension of the array and the precision to work on. If one of these arguments is not present, the program will detect and halt. Example:
mpirun -np 5 ./executable 16 0.04. The first argument after the executable program name will be the dimension and the second will the precision.

NOTE: This is a PROOF OF CONCEPT for demeostration, not an absolute state of the art implementation of the problem.
