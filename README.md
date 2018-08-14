This is an implementation of the decremental dynamic algorithm to compute mutually connected components
by S. Hwang, S. Choi, D. Lee and B. Kahng (Physical Review E, 91, 022814 (2015)).

`Implementation` folder contains the implementation as a header only library.
`Example` folder contains some examples using it.
`Test` folder contains tests for the library.

Reading `Example/DecrementalGMCC.cpp` and `Test/DecrementalMCCTest.cpp` may help you to get how to use it.
The implementation uses the C++11 standard, so you may need to set a compiler flag  (e.g. `-std=c++11` or `-std=c++1z` for gcc or clang) to use it.
It has a compile time dependency on Boost.Optional.