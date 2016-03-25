#pragma once

#include <papi.h>
#include <omp.h>
#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <fstream>

namespace pompi
{

    enum OutputFormat
    {
        GNUPLOT,
    };


    class Base
    {

        private:

            double execution_time_start_;
            double execution_time_end_;
            std::vector< int > papi_events_;
            std::vector< std::string > papi_event_names_;
            std::vector< std::vector< long long > > thread_data_;
            int max_threads_;

        public:

            Base();

            void AddEvent(char * event);

            void Start();

            void Stop();

            double GetExecutionTime();

            double GetAverageExecutionTime(int trials);

            void ClearAllCounters();

            void ClearTimers();

            void PrintThreadResultsToFile(int thread_id, int total_threads, char * file_name, OutputFormat format);

            void PrintThreadResults(int thread_id, int total_threads);

            void PrintAggregatedResultsToFile(int total_threads, char * file_name, OutputFormat format);

            void PrintAggregatedResults(int total_threads);

        private:

            void GetThreadCounters(int thread_id, long long * counters);

            void GetAggregatedCounters(long long * counters);

            void ClearThreadCounters(int thread_id);

            void PrintThreadGnuplot(int thread_id, int total_threads, std::ofstream &output);

            void PrintAggregatedGnuplot(int total_threads, std::ofstream &output);
    };



    Base::Base()
    {
        max_threads_ = omp_get_max_threads();
        std::clog << "[Log] Maximum number of threads is " << max_threads_ << std::endl;

        thread_data_.resize(max_threads_);
        ClearTimers();

        int PAPI_return_value = PAPI_library_init(PAPI_VER_CURRENT);
        if(PAPI_return_value != PAPI_VER_CURRENT)
        {
            std::cerr << "[Error] Could not initialize PAPI library" << std::endl;
            exit(1);
        }
        std::clog << "[Log] PAPI library successfully initialized" << std::endl;

        PAPI_return_value = PAPI_thread_init((unsigned long (*)(void)) (omp_get_thread_num));
        if(PAPI_return_value != PAPI_OK)
        {
            std::cerr << "[Error] Coult not initialize OpenMP threading for PAPI" << std::endl;
            exit(1);
        }
        std::clog << "[Log] OpenMP threading for PAPI successfully initialized" << std::endl;

        // Resolving PAPI events from PAPI_EVENTS enviroment variable
        char * papi_events = getenv("PAPI_EVENTS");
        char * event = strtok(papi_events, "|");
        int event_id;
        while(event != NULL)
        {
            AddEvent(event);
            event = strtok(NULL, "|");
        }
    }

    void Base::AddEvent(char * event)
    {
        int event_id;
        int PAPI_return_value = PAPI_event_name_to_code(event, &event_id);
        if(PAPI_return_value != PAPI_OK)
            std::cerr << "[Warning] Papi event `" << event << "` does not exists, skipping" << std::endl;
        else
        {
            if(std::find(papi_events_.begin(), papi_events_.end(), event_id) == papi_events_.end())
            {
              papi_event_names_.push_back(std::string(event));
              papi_events_.push_back(event_id);
              for(int i = 0; i < max_threads_; ++i)
                thread_data_[i].push_back(0);
              std::clog << "[Log] Adding papi event `" << event << "`" << std::endl;
            }
            else
                std::cerr << "[Warning] Papi event `" << event << "` already listed, skipping" << std::endl;
        }
    }

    void Base::Start()
    {
        int thread_id = omp_get_thread_num();

        int PAPI_return_value = PAPI_start_counters(&papi_events_[0], papi_events_.size());
        if(PAPI_return_value != PAPI_OK)
        {
            #pragma omp critical
            std::cerr << "[Error] Could not start counters in thread #" << thread_id << std::endl;
            exit(1);
        }

        #pragma omp single
        {
            if(execution_time_start_ == 0)
                execution_time_start_ = omp_get_wtime();
        }
    }

    void Base::Stop()
    {
        #pragma omp single
        execution_time_end_ = omp_get_wtime();

        long long counter_values[papi_events_.size()];

        int thread_id = omp_get_thread_num();
        int PAPI_return_value = PAPI_stop_counters(counter_values, papi_events_.size());
        if(PAPI_return_value != PAPI_OK)
        {
            #pragma omp critical
            std::cerr << "[Error] Could not stop counters in thread #" << thread_id << std::endl;
            exit(1);
        }

        for(int event = 0; event < papi_events_.size(); ++event)
            thread_data_[thread_id][event] += counter_values[event];

        PAPI_unregister_thread();
    }

    double Base::GetExecutionTime()
    {
        return execution_time_end_ - execution_time_start_;
    }

    double Base::GetAverageExecutionTime(int trials)
    {
        return GetExecutionTime()/trials;
    }

    void Base::ClearAllCounters()
    {
        for(int thread = 0; thread < max_threads_; ++thread)
            ClearThreadCounters(thread);

    }

    void Base::ClearTimers()
    {
        execution_time_start_ = 0;
        execution_time_end_ = 0;
    }

    void Base::GetThreadCounters(int thread_id, long long * counters)
    {
        for(int event = 0; event < papi_events_.size(); ++event)
            counters[event] = thread_data_[thread_id][event];
    }

    void Base::GetAggregatedCounters(long long * counters)
    {
        long long thread_counters[papi_events_.size()];
        for(int event = 0; event < papi_events_.size(); ++event)
            counters[event] = 0;

        for(int thread = 0; thread < max_threads_; ++thread)
        {
            GetThreadCounters(thread, thread_counters);
            for(int event = 0; event < papi_events_.size(); ++event)
                counters[event] += thread_counters[event];
        }
    }

    void Base::ClearThreadCounters(int thread_id)
    {
        for(int event = 0; event < papi_events_.size(); ++event)
            thread_data_[thread_id][event] = 0;
    }

    void Base::PrintThreadResults(int thread_id, int total_threads)
    {
        if((thread_id > max_threads_) || (thread_id < 0))
        {
            std::cerr << "[Warning] invalid thread_id in PrintThreadResults" << std::endl;
            return;
        }

        long long results[papi_events_.size()];
        GetThreadCounters(thread_id, results);

        std::cout << "Results for thread #" << thread_id << " out of " << total_threads << " threads" << std::endl;
        for(int event = 0; event < papi_events_.size(); ++event)
        {
            std::cout << papi_event_names_[event] + ':'
                      << std::setw(20) << results[event]
                      << std::endl;
        }
    }

    void Base::PrintAggregatedResults(int total_threads)
    {
        std::cout << "Aggregated results on " << total_threads << " threads" << std::endl;
        long long results[papi_events_.size()];
        GetAggregatedCounters(results);

        for(int event = 0; event < papi_events_.size(); ++event)
        {
            std::cout << papi_event_names_[event] + ':'
                      << std::setw(20) << results[event]
                      << std::endl;
        }
    }

    void Base::PrintThreadResultsToFile(int thread_id, int total_threads, char * file_name, OutputFormat format)
    {
        if((thread_id < 0)||(thread_id > total_threads))
        {
            std::cerr << "[Warning] Unable to print thread result to file" << std::endl;
            return;
        }

        std::ofstream output;
        output.open(file_name, std::ofstream::app);

        switch(format)
        {
            case GNUPLOT: {
                PrintThreadGnuplot(thread_id, total_threads, output);
                break;
            }
        }

        output.close();
    }

    void Base::PrintAggregatedResultsToFile(int total_threads, char * file_name, OutputFormat format)
    {
        std::ofstream output;
        output.open(file_name, std::ofstream::app);

        switch(format)
        {
            case GNUPLOT: {
                PrintAggregatedGnuplot(total_threads, output);
                break;
            }

        }

        output.close();
    }

    void Base::PrintThreadGnuplot(int thread_id, int total_threads, std::ofstream &output)
    {
        output << std::setw(8) << "#THREADS";
        for(int event = 0; event < papi_events_.size(); ++event)
            output << std::setw(16) << papi_event_names_[event];
        output << std::endl;

        long long results[papi_events_.size()];
        GetThreadCounters(thread_id, results);

        output << std::setw(8) << total_threads;
        for(int event = 0; event < papi_events_.size(); ++event)
            output << std::setw(16) << results[event];
        output << std::endl;
    }

    void Base::PrintAggregatedGnuplot(int total_threads, std::ofstream &output)
    {
        output << std::setw(8) << "#THREADS";
        for(int event = 0; event < papi_events_.size(); ++event)
            output << std::setw(16) << papi_event_names_[event];
        output << std::endl;

        long long results[papi_events_.size()];
        GetAggregatedCounters(results);

        output << std::setw(8) << total_threads;
        for(int event = 0; event < papi_events_.size(); ++event)
            output << std::setw(16) << results[event];
        output << std::endl;
    }
}
