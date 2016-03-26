POMPI - PAPI for OpenMP
=======================

POMPI provides all necessary tools for successful monitoring of OpenMP multithreaded applications using Performance Application Programmable Interface (PAPI) in C++ code. Motivation for creation of this rather simple header file was the need for a working PAPI wrapper in multithreaded enviroment, that would be easy to use, well documented and open source.

Requirements
============

* Working PAPI. Check out [official PAPI documentation page](http://icl.cs.utk.edu/projects/papi/wiki/Main_Page). Compile with `-lpapi`
* Wrapper uses OpenMP runtime functions so don't forget to compile with `-fopenmp`. If you want to learn more about OpenMP and parallelism check out [OpenMP pages](http://openmp.org/wp/)

Usage
=====

To make use of POMPI just copy `pompi.h` into your project directory and include it. This will provide you with `pompi::` namespace. Initialization of POMPI is done by instantiating `pompi::Base` class. All methods are then accessed by this instance. See documentation
for all available methods.

The most important are:

* `AddEvent(char *)` - adds new PAPI event to be monitored.
* `Start()` - starts performance monitoring.
* `Stop()` - stop performance monitoring.

There are several rules to using POMPI:

* There should be only one instance of `pompi::Base` at a time (singleton rework comming soon).
* Both `Start()` and `Stop()` methods must be called inside of the same parallel region.
* A call to `Start()` must precede call to `Stop()`.

Example usage
-------------

``` c++
pompi::Base pompi_instance;
pompi_instance.AddEvent("PAPI_L2_TCA");
pompi_instance.AddEvent("PAPI_L2_TCM");

#pragma omp parallel num_threads(x)
{
  pompi_instance.Start();
  insanely_interesting_cache_operation();
  pompi_instance.Stop();
}

std::cout << "Execution time: " << pompi_instance.GetExecutionTime() << std::endl;
pompi_instance.PrintResults(x);
```

Benchmarking
------------

POMPI can be used for benchmarking as it is able to compute derived events such as L2 miss ratio. It can provide results for each separate thread as well as aggregated results on all threads. Timers in POMPI are made with support for benchmarking as well. POMPI is also parsing enviroment variable `PAPI_EVENTS` and adding these in the constructor. Lets see modified example:

``` c++
pompi::Base pompi_instance;
pompi_instance.AddEvent("PAPI_L2_TCA");
pompi_instance.AddEvent("PAPI_L2_TCM");

for(int i = 0; i < REPEAT_COUNT; ++i)
{
  #pragma omp parallel num_threads(x)
  {
    pompi_instance.Start();
    insanely_interesting_cache_operation();
    pompi_instance.Stop();
  }
}
```

In this example, we are repeating same operation over and over again to gain better precision. POMPI, instead of saving data from the last iteration of the loop, is accumulating data as well as execution time. In the end, we can use `GetAverageExecutionTime(int REPEAT_COUNT)` to get average time. Results of `Print` methods are accumulated as well.

Testing
=======

POMPI has been tested using Intel C compiler on up to 24 threads successfully. Feel free to test it further.

Documentation
=============


	git clone https://github.com/LukasSlouka/pompi.git
	cd pompi
	doxygen Doxyfile
	firefox html/index.html (or something)


Further development
===================

For now there are only 3 cache oriented derived events available as well as only one file output format supported (GNUPLOT). Adding more possibilities should not be hard at all so feel free to contribute as you need it.

