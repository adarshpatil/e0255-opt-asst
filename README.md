E0255 Assignment

Deadline: Apr 26, 2015 11:59pm

The objective of the assignment is to optimize the Harris corner 
detection algorithm using transformations learnt during the course, in 
particular, locality optimizations, vectorization, and parallelization.

Put your function in a separate file harris.opt.cpp and name it 
harris_opt.  With the video demo, hitting the key 'o' should switch from 
the base implementation provided in harris.cpp (harris_base) to yours 
(harris_opt).  You'll have to update video_benchmark.py slightly to 
allow this.

All performance measurements should be taken on CL1 workstations. The 
configuration of the workstations is 

HP Z230 Tower Workstation
Intel(R) Core(TM) i7-4770 CPU @ 3.40GHz - (Haswell microarchitecture)
4 physical cores  (x 2 hyperthreads / core)
64 KB L1 private / 256 KB L2 private / 8 MB L3 shared cache
32 GB DDR3-1600 memory
2 x 1 TB storage in RAID 1
128GB SSD storage (for root file system)
NVIDIA Quadro K620 (384 cores, 2GB GDDR3 RAM, PCI-express card) 


WHAT TO SUBMIT
---------------
Submit (1) your harris.opt.cpp with the modified function named as 
described above, and with the same signature as harris_base, (2) the 
modified video_benchmark.py, and (3) a file named 'report.txt' with a 
description of the optimizations you performed, an explanation of why 
you performed those optimizations, and the performance you measured 
while running on 1, 2, 3, and 4 cores (with no hyperthreads in any case) 
on a CL1 workstation. All three files should be sent in a single tar 
gzipped file named <your name>.tar.gz. 


Compiling the reference implementation
--------------------------------------
Install OpenCV with QT or gtk

These two urls should help. 

http://karytech.blogspot.in/2012/05/opencv-24-on-ubuntu-1204.html

http://opencv-python-tutroals.readthedocs.org/en/latest/py_tutorials/py_setup/py_setup_in_fedora/py_setup_in_fedora.html

One can also install from source using the documentation on the OpenCV website

icpc -xhost -openmp -O3 harris.cpp -L /usr/local/lib/ -lopencv_imgproc -lopencv_core -lopencv_highgui -o harris -DANALYZE -DSHOQ -DNRUNS=5

g++ -openmp -O3 harris.cpp -L /usr/local/lib/ -lopencv_imgproc -lopencv_core -lopencv_highgui -o harris -DANALYZE -DSHOW -DNRUNS=5

./harris.out path_to_image_file

Running the python demo
-----------------------
Install Python, NumPy (the urls above cover this installation) 

Compile the reference implementation in the video folder as a shared library

icpc -xhost -openmp -fPIC -shared -o harris.so harris_extern.cpp

g++ -openmp -fPIC -shared -o harris.so harris_extern.cpp

Run the python script as

python video_benchmark.py path_to_video_file 1/0

1 - for displaying the video

0 - for running 25 frames using OpenCV and 25 frames using reference implementation

Once the video is running
 
h - toggles harris mode on/off

space - toggles between opencv and base or optimized implementation

o - toggles between base and optimized implementations

## Optimizations applied
### Compiler flags
* -O3 
* -fprefetch-loop-arrays
* -fopenmp
* -ffast-math
* -fprofile-generate
* -fprofile-use

### Code optimization
* Tiling 
* Tile sized intermediate Scratchpad arrays 
* Unroll Jam i loop

### Pragma Used
* pragma omp parallel for
* pragma GCC ivdep

### Evaluation
#### Image size
* 21600px X 10800px
* 29MB
* Best reference time vs Worst optimized time measured

#### Processor
* Intel(R) Core(TM) i7-4770 CPU @ 3.40GHz
* intel_pstate driver; performance governer
* Single socket, 4 physical cores, hyperthreading disabled
* Caches 
    * L1d cache:             32K
    * L1i cache:             32K
    * L2 cache:              256K
    * L3 cache:              8192K

#### Performance Speedup
GCC 4.9.2
* Single Core + vectorize = **3.29x**
* Multi Core + vectorize = **11.34x**
Detailed performance speedup comparison of ICC vs GCC and vectorization, parallelism etc. available in report.pdf

