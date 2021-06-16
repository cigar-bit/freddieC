# freddieC
Rough C++ implementation of freddie_segment

To build install cmake if not already installed and run the following:

$ cmake -DCMAKE_BUILD_TYPE=Release .

$ make

Note: If you prefer building in another directory that is fine as well just move the final freddie_segment executable into the parent directory of the test/ directory

To run:

$ ./freddie_segment -s test/split/ --outdir test/freddie_segment/ -t 1

Arguments are the same as the original Freddie, this implementation can support multithreading as well.
