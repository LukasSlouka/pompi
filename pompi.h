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

    #define ALL_THREADS -1

    enum OutputFormat
    {
        GNUPLOT,
    };


    enum PapiDerivedStat
    {
        D_L1_TMR,
        D_L2_TMR,
        D_L3_TMR,
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

            void ClearCounters(int thread_id = ALL_THREADS);

            void ClearTimers();

            void PrintResults(int total_threads, int thread_id = ALL_THREADS);

            void PrintResultsToFile(int total_threads, char * file_name, OutputFormat format, int thread_id = ALL_THREADS);

        private:

            void GetCounters(long long * counters, int thread_id = ALL_THREADS);

            int GetEventIndex(int event_code);

            bool EventAvailable(int event_code);

            void GetDerivedStats(std::vector< PapiDerivedStat > & stats);

            std::string GetDerivedStatName(PapiDerivedStat stat);

            void PrintGnuplot(int total_threads, std::ofstream &output, int thread_id = ALL_THREADS);

            double ComputeDerivedStat(PapiDerivedStat stat, int thread_id = ALL_THREADS);
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


    void Base::ClearCounters(int thread_id)
    {
        if((thread_id < 0)||(thread_id > max_threads_))
        {
            for(int event = 0; event < papi_events_.size(); ++event)
                thread_data_[thread_id][event] = 0;
        }
        else
        {
            for(int thread = 0; thread < max_threads_; ++thread)
            for(int event = 0; event < papi_events_.size(); ++event)
                thread_data_[thread][event] = 0;                        
        }
    }


    void Base::ClearTimers()
    {
        execution_time_start_ = 0;
        execution_time_end_ = 0;
    }


    void Base::GetCounters(long long * counters, int thread_id)
    {
        for(int event = 0; event < papi_events_.size(); ++event)
            counters[event] = 0;

        if(thread_id == ALL_THREADS)
        {
            for(int thread = 0; thread < max_threads_; ++thread)
            for(int event = 0; event < papi_events_.size(); ++event)
                counters[event] += thread_data_[thread][event];
        }
        else
        {
            for(int event = 0; event < papi_events_.size(); ++event)
                counters[event] = thread_data_[thread_id][event];
        }
        
    }


    bool Base::EventAvailable(int event_code)
    {
        std::vector<int>::iterator event;
        event = std::find(papi_events_.begin(), papi_events_.end(), event_code);
        return (event == papi_events_.end()) ? false : true;
    }


    int Base::GetEventIndex(int event_code)
    {
        for(int index = 0; index < papi_events_.size(); ++index)
            if(papi_events_[index] == event_code)
                return index;
        return 0;
    }


    void Base::GetDerivedStats(std::vector< PapiDerivedStat > &stats)
    {
        if (EventAvailable(PAPI_LD_INS) && EventAvailable(PAPI_SR_INS) && EventAvailable(PAPI_L1_TCM))
            stats.push_back(D_L1_TMR);
        
        if (EventAvailable(PAPI_L2_TCA) && EventAvailable(PAPI_L2_TCM))
            stats.push_back(D_L2_TMR);
        
        if (EventAvailable(PAPI_L3_TCA) && EventAvailable(PAPI_L3_TCM))
            stats.push_back(D_L3_TMR);
    }


    double Base::ComputeDerivedStat(PapiDerivedStat stat, int thread_id)
    {
        long long counters[papi_events_.size()];
        GetCounters(counters, thread_id);

        switch(stat)
        {
            case D_L1_TMR: {
                return (counters[GetEventIndex(PAPI_L1_TCM)] / (double)(counters[GetEventIndex(PAPI_SR_INS)] + counters[GetEventIndex(PAPI_LD_INS)]));
            }
            case D_L2_TMR: {
                return (counters[GetEventIndex(PAPI_L2_TCM)] / (double)counters[GetEventIndex(PAPI_L2_TCA)]);
            }
            case D_L3_TMR: {
                return (counters[GetEventIndex(PAPI_L3_TCM)] / (double)counters[GetEventIndex(PAPI_L3_TCA)]);
            }
        }
    }


    std::string Base::GetDerivedStatName(PapiDerivedStat stat)
    {
        switch(stat)
        {
            case D_L1_TMR: {
                return "D_L1_TMR";
            }
            case D_L2_TMR: {
                return "D_L2_TMR";
            }
            case D_L3_TMR: {
                return "D_L3_TMR";
            }
        }
    }


    void Base::PrintResults(int total_threads, int thread_id)
    {
        std::vector< PapiDerivedStat > stats;
        GetDerivedStats(stats);

        if((thread_id >= 0)&&(thread_id <= max_threads_))
            std::cout << "Results for thread #" << thread_id << " out of "
                      << total_threads << " threads" << std::endl;
        else
        {
            std::cout << "Aggregated results on " << total_threads << " threads" << std::endl;
            thread_id = ALL_THREADS;
        }

        long long results[papi_events_.size()];
        GetCounters(results, thread_id);

        for(int event = 0; event < papi_events_.size(); ++event)
        {
            std::cout << papi_event_names_[event] + ':'
                      << std::setw(20) << results[event]
                      << std::endl;
        }

        for(int d_event = 0; d_event < stats.size(); ++d_event)
        {
            std::cout << GetDerivedStatName(stats[d_event]) + ':'
                      << std::setw(20) << ComputeDerivedStat(stats[d_event], thread_id)
                      << std::endl;
        }
    }


    void Base::PrintResultsToFile(int total_threads, char * file_name, OutputFormat format, int thread_id)
    {
        if((thread_id < 0)||(thread_id > total_threads))
            thread_id = ALL_THREADS;

        std::ofstream output;
        output.open(file_name, std::ofstream::app);

        switch(format)
        {
            case GNUPLOT: {
                PrintGnuplot(total_threads, output, thread_id);
                break;
            }
        }

        output.close();
    }


    void Base::PrintGnuplot(int total_threads, std::ofstream &output, int thread_id)
    {
        std::vector< PapiDerivedStat > stats;
        GetDerivedStats(stats);

        if(thread_id == ALL_THREADS)
            output << std::setw(8) << "#THREADS";
        else
            output << std::setw(7) << "#THREAD";

        for(int event = 0; event < papi_events_.size(); ++event)
            output << std::setw(16) << papi_event_names_[event];
        for(int d_event = 0; d_event < stats.size(); ++d_event)
            output << std::setw(16) << GetDerivedStatName(stats[d_event]);
        output << std::endl;

        long long results[papi_events_.size()];
        GetCounters(results, thread_id);

        if(thread_id == ALL_THREADS)
            output << std::setw(8) << total_threads;
        else
            output << std::setw(7) << thread_id;

        for(int event = 0; event < papi_events_.size(); ++event)
            output << std::setw(16) << results[event];
        for(int d_event = 0; d_event < stats.size(); ++d_event)
            output << std::setw(16) << ComputeDerivedStat(stats[d_event], thread_id);
        output << std::endl;
    }
}
