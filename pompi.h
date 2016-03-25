#pragma once

#include <papi.h>
#include <omp.h>
#include <vector>
#include <string>
#include <iostream>


namespace pompi
{

    // Declaration

    class Base
    {

        private:

            double execution_time_start_;
            double execution_time_end_;
            std::vector< int > papi_events_;
            std::vector< std::vector< long long > > thread_data_;
            int max_threads_;


        public:
            Base();
    };

    // Implementation

    Base::Base()
    {
        max_threads_ = omp_get_max_threads();
        std::clog << "Maximum number of threads is: " << max_threads_ << std::endl;

        thread_data_.resize(max_threads_);

        int PAPI_return_value;

        PAPI_return_value = PAPI_library_init(PAPI_VER_CURRENT);
        if(PAPI_return_value != PAPI_VER_CURRENT)
        {
            std::cerr << "Could not initialize PAPI library" << std::endl;
            exit(1);
        }
        std::clog << "PAPI library successfully initialized" << std::endl;

        PAPI_return_value = PAPI_thread_init((unsigned long (*)(void)) (omp_get_thread_num));
        if (PAPI_return_value != PAPI_OK)
        {
            std::cerr << "Coult not initialize OpenMP threading for PAPI" << std::endl;
            exit(1);
        }
        std::clog << "OpenMP threading for PAPI successfully initialized" << std::endl;
    }

}
