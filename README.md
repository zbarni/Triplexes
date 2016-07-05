# triplexator

1. Compiling
    
    1.1 binary: 

        *) mkdir -p build/Release && cd build/Release
        *) cmake ../.. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -G "Unix Makefiles"
        *) make

    1.2 shared library (for python bindings):
    
        *) mkdir -p build/Release && cd build/Release
        *) cmake ../.. -DSHAREDLIBRARY=TRUE -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -G "Unix Makefiles"
        *) make

2. Python bindings

    *) set your TRIPLEXATOR_LIBRARY environment variable to the location of the Triplexator shared library (libtriplexator.so)
    
    *) call runTriplexator function with the regular Triplexator parameters as string param.

3. Run Triplexator with the following options
    
    - **--bit-parallel**: for bit parallel (faster) version of brute-force
    - **--bit-parallel-local**: for semi-palindrom search
