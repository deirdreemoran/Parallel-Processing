# Parallel-Programming
Parallel Programming Projects 

Simple openMP:  experiment that takes large number arrays and multiples them together 
in order to test parallel processing execution times.

Open MP-mem-alloc-grain:  OpenMP N-body Problem - Examines variations in performance based on static/dynamic and fine/coarse grained parallelism.

OpenMP-Numeric-integration:  Numeric Integration with OpenMP.  Calculates the volume between two Bézier surfaces 
and tests the parallel processing execution times based on size of set of 
subdivisions and number of threads used.

FalseSharing:  In False Sharing, two (or more) cores are writing to variables that live in the same cache line. 
As soon as each one writes to the cache line, it invalidates the contents of the other, 
requiring a time-consuming re-load.  For this project, we avoid false-sharing by implementing
two fixes.  The first fix, project3_fix1.cpp, uses contiguous array padding to ensure that
each individual block starts and ends on a cache boundary.  The second fix, project3_fix2.cpp, uses
a temporary private variable (created in each core’s stack) so that very little to no cache line 
conflicts occur.
