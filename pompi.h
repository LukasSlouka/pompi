/*
    The MIT License (MIT)

    Copyright (c) 2016 Lukáš Slouka

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/


/**
 * @file pompi.h
 * @brief This file contains all necessary tools for successful monitoring
 *        of OpenMP multithreaded applications using Performance Application
 *        Programmable Interface (PAPI) in C++ code. For all intents and purposes
 *        it acts like a wrapper.
 * @author Lukas Slouka <lukslouka@gmail.com>
 * @version 1.0
 */

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

/// pompi namespace.
namespace pompi
{

    /// Used in methods that support output on either only one thread or all threads.
    #define ALL_THREADS -1

    // Add new output formats here
    /**
     * Output formats enum.
     * Contains all available kinds of output format.
     */
    enum OutputFormat
    {
        GNUPLOT,    /**< Formats performance data to be usable with GNUplot. */
    };

    // Add new derived stats here
    /**
     * Perived statistics available from PAPI counters.
     * Contains all available statistics derived from standard PAPI events.
     */
    enum PapiDerivedStat
    {
        D_L1_TMR,   /**< Derived L1 cache total missrate [%]*/
        D_L2_TMR,   /**< Derived L2 cache total missrate [%]*/
        D_L3_TMR,   /**< Derived L3 cache total missrate [%]*/
        D_L1_DMR,   /**< Derived L1 cache data missrate [%]*/
        D_L2_DMR,   /**< Derived L2 cache data missrate [%]*/
        D_L3_DMR,   /**< Derived L3 cache data missrate [%]*/
        D_BR_MPR,   /**< Derived branch missprediction rate [%]*/
    };

    /**
     * Pompi Base class.
     * This class is heart and soul of pompi wrapper. Contains all methods
     * and attributes needed for successful monitoring of performance.
     */
    class Base
    {

        private:

            /**
             * Execution time start variable.
             * Contains time of execution begining.
             * @warning Note that it does not change value with each call to Start method. 
             *          Its stays the same until it is cleared by ClearTimers method.
             * This behaviour is due to use of repeating the same action multiple
             * times in benchmarking, in which case its needed to accumulate execution time
             * rather than save only the execution time of the last iteration.
             * @see ClearTimers()
             * @see Start()
             */
            double execution_time_start_;

            /**
             * Execution time stop variable.
             * Stop time is updated with each call to Stop method.
             * @see Stop()
             * @see ClearTimers()
             */
            double execution_time_end_;

            /**
             * Vector of papi event codes.
             * Dependant on PAPI_EVENTS enviroment variable.
             */
            std::vector< int > papi_events_;

            /**
             * Vector of papi event names.
             * Content corresponds with papi_events_.
             */
            std::vector< std::string > papi_event_names_;

            /**
             * Vector of thread hardware counter values.
             * Filled by each thread after call to Stop.
             * @see ClearCounters()
             * @see GetCounters()
             * @see Stop()
             */
            std::vector< std::vector< long long > > thread_data_;

            /**
             * Maximum number of threads.
             * Initialized by call to omp_get_max_threads().
             */
            int max_threads_;

        public:

            /**
             * Base constructor
             * Initializes PAPI library and PAPI library threading
             * support for OpenMP. Parses PAPI_EVENTS enviroment variable
             * and modifies papi_events_ and papi_event_names_ accordingly.
             * @warning Should be instantiated only once.
             */
            Base();

            /**
             * Tells pompi to monitor additional event.
             * @param event C string name of PAPI event
             */
            void AddEvent(char * event);

            /**
             * Starts counting of all monitored events on all threads.
             * This method is complemented by Stop method.
             * @warning Must be called inside of a OpenMP parallel region
             * @see Stop()
             */
            void Start();

            /**
             * Stops counting of all monitored events on all threads
             * This method is complemented by Star method.
             * @warning Must be called inside of a OpenMP parallel region
             * @see Start()
             */
            void Stop();

            /**
             * Provides execution time between first call to Start() and
             * last call to Stop(). This behaviour is reset by ClearTimers().
             * @return execution time
             * @see Start()
             * @see Stop()
             * @see ClearTimers()
             */
            double GetExecutionTime();

            /**
             * Similar to GetExecutionTime. Used to calculate execution time 
             * averaged on number of trials.
             * @param  trials repetition count.
             * @return        average execution time.
             * @see GetExecutionTime()
             */
            double GetAverageExecutionTime(int trials);

            /**
             * Sets all counters to zero. Number of counter depends on
             * maximum possible threads and count of monitored papi events.
             * Can be used to clear counters of any singular thread or all of them.
             * @param thread_id Integer ID of a thread. If thread with thread_id
             *        does not exist, all counters are cleared.
             */
            void ClearCounters(int thread_id = ALL_THREADS);

            /**
             * Sets both execution time timers to zero.
             */
            void ClearTimers();

            /**
             * Prints results on standard output. Can be used to output counters
             * of any singular thread or all threads aggregated (implicit).
             * @param total_threads Total number of threads used during execution.
             * @param thread_id     ID of a thread to be outputted (invalid thread ID
             *                      results in aggregated output on all threads).
             */
            void PrintResults(int total_threads, int thread_id = ALL_THREADS);

            /**
             * Simillar to PrintResults. Takes additional file description 
             * parameters to define output file.
             * @param total_threads Total number of threads used during execution.
             * @param file_name     C string name of a output file.
             * @param format        Output format.
             * @param thread_id     ID of a thread to be outputted (invalid thread ID
             *                      results in aggregated output on all threads).
             */
            void PrintResultsToFile(int total_threads, char * file_name, OutputFormat format, int thread_id = ALL_THREADS);

        private:

            /**
             * Private method used to obtain counters
             * @param counters  array of values where counters are to be saved.
             * @param thread_id Counters of a thread with this ID will be extracted.
             */
            void GetCounters(long long * counters, int thread_id = ALL_THREADS);

            /**
             * Private method used to get index of an event within papi_events_.
             * @param  event_code Integer event code.
             * @return            index.
             */
            int GetEventIndex(int event_code);

            /**
             * Private predicate checking availablity of PAPI event within
             * monitored events.
             * @param  event_code Integer event code.
             * @return            Boolean result of predicate.
             */
            bool EventAvailable(int event_code);

            /**
             * Private method obtaining all derived event that are available.
             * @param stats Vector filled with derived stats.
             */
            void GetDerivedStats(std::vector< PapiDerivedStat > & stats);

            /**
             * Private method translating derived stat code to string.
             * @param  stat Derived stat code.
             * @return      String representation of derived stat.
             */
            std::string GetDerivedStatName(PapiDerivedStat stat);

            /**
             * Private method that computes derived stat value.
             * @param  stat      Derived stat to be computed
             * @param  thread_id Derived stat will be computed for a thread with
             *                   this thread_id. In case thread_id is invalid,
             *                   the stat will be computed aggregated for all threads.
             * @return           Derived stat value
             */
            double ComputeDerivedStat(PapiDerivedStat stat, int thread_id = ALL_THREADS);

            /**
             * Private method defining print format for GNUPLOT
             * @param total_threads Total number of threads used during execution
             * @param output        Output stream
             * @param thread_id     Output will be performed for thread with this ID.
             *                      In case thread_id is invalid, output will be
             *                      performed aggregated for all threads.
             */
            void PrintGnuplot(int total_threads, std::ofstream &output, int thread_id = ALL_THREADS);
    };

    
    /////////////////////////////////////
    //  Implementation PUBLIC METHODS  //
    /////////////////////////////////////

    
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
            for(int thread = 0; thread < max_threads_; ++thread)
            for(int event = 0; event < papi_events_.size(); ++event)
                thread_data_[thread][event] = 0;   
        }
        else
        {
            for(int event = 0; event < papi_events_.size(); ++event)
                thread_data_[thread_id][event] = 0;                     
        }
    }


    
    void Base::ClearTimers()
    {
        execution_time_start_ = 0;
        execution_time_end_ = 0;
    }

    //////////////////////////////////////
    //  Implementation PRIVATE methods  //
    //////////////////////////////////////

    
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


    // Add new derived stats here
    void Base::GetDerivedStats(std::vector< PapiDerivedStat > &stats)
    {
        if (EventAvailable(PAPI_LD_INS) && EventAvailable(PAPI_SR_INS) && EventAvailable(PAPI_L1_TCM))
            stats.push_back(D_L1_TMR);
        
        if (EventAvailable(PAPI_L2_TCA) && EventAvailable(PAPI_L2_TCM))
            stats.push_back(D_L2_TMR);
        
        if (EventAvailable(PAPI_L3_TCA) && EventAvailable(PAPI_L3_TCM))
            stats.push_back(D_L3_TMR);

        if (EventAvailable(PAPI_L1_DCA) && EventAvailable(PAPI_L1_DCM))
            stats.push_back(D_L1_DMR);

        if (EventAvailable(PAPI_L2_DCA) && EventAvailable(PAPI_L2_DCM))
            stats.push_back(D_L2_DMR);

        if (EventAvailable(PAPI_L3_DCA) && EventAvailable(PAPI_L3_DCM))
            stats.push_back(D_L3_DMR);

        if (( EventAvailable(PAPI_BR_MSP) && EventAvailable(PAPI_BR_CN) )||
            ( EventAvailable(PAPI_BR_MSP) && EventAvailable(PAPI_BR_PRC) ))
            stats.push_back(D_BR_MPR);
    }


    double Base::ComputeDerivedStat(PapiDerivedStat stat, int thread_id)
    {
        long long counters[papi_events_.size()];
        GetCounters(counters, thread_id);

        // Add new derived stat computations here
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
            case D_L1_DMR: {
                return (counters[GetEventIndex(PAPI_L1_DCM)] / (double)counters[GetEventIndex(PAPI_L1_DCA)]);
            }
            case D_L2_DMR: {
                return (counters[GetEventIndex(PAPI_L2_DCM)] / (double)counters[GetEventIndex(PAPI_L2_DCA)]);
            }
            case D_L3_DMR: {
                return (counters[GetEventIndex(PAPI_L3_DCM)] / (double)counters[GetEventIndex(PAPI_L3_DCA)]);
            }
            case D_BR_MPR: {
                if(EventAvailable(PAPI_BR_CN))
                    return (counters[GetEventIndex(PAPI_BR_MSP)] / (double)counters[GetEventIndex(PAPI_BR_CN)]);
                else
                    return (counters[GetEventIndex(PAPI_BR_MSP)] / (double)(counters[GetEventIndex(PAPI_BR_MSP)] + counters[GetEventIndex(PAPI_BR_PRC)]));
            }
        }
    }

    // Add new derived stat translations here
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
            case D_L1_DMR: {
                return "D_L1_DMR";
            }
            case D_L2_DMR: {
                return "D_L2_DMR";
            }
            case D_L3_DMR: {
                return "D_L3_DMR";
            }
            case D_BR_MPR: {
                return "D_BR_MPR";
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

        // Add new output format methods here
        switch(format)
        {
            case GNUPLOT: {
                PrintGnuplot(total_threads, output, thread_id);
                break;
            }
        }

        output.close();
    }


    // Format specific print methods

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
