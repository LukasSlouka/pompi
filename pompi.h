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
     * Pompi timer class
     * Timer class that can be used by itself, but is part of wrapper
     */
    class PompiTimer
    {
        private:

            /**
             * First call to timer
             */
            double first_call_time_;

            /**
             * Last call to timer
             */
            double last_call_time_;

            /**
             * Inner timer call (does not change first call)
             */
            double inner_call_time_;

            /**
             * Inner calls aggregated duration
             */
            double inner_calls_sum_;

            /**
             * Inner calls count
             */
            unsigned int repetition_count_;


        public:

            /**
             * A constructor
             * Initializes all counters to zero
             */
            PompiTimer()
            {
                first_call_time_ = 0.0f;
                last_call_time_ = 0.0f;
                inner_call_time_ = 0.0f;
                inner_calls_sum_ = 0.0f;
                repetition_count_ = 0.0f;
            }

            /**
             * Begins timing
             * Sets first call if it is zero, otherwise sets inner call time
             */
            void BeginTiming()
            {
                inner_call_time_ = omp_get_wtime();
                if(first_call_time_ == 0.0f)
                    first_call_time_ = inner_call_time_;
            }

            /**
             * Ends timing
             */
            int EndTiming()
            {
                last_call_time_ = omp_get_wtime();
                repetition_count_++;
                inner_calls_sum_ += last_call_time_ - inner_call_time_;
            }

            /**
             * Clears entire timer
             */
            void inline ResetTimer()
            {
                first_call_time_ = 0.0f;
                last_call_time_ = 0.0f;
                inner_call_time_ = 0.0f;
                inner_calls_sum_ = 0.0f;
                repetition_count_ = 0.0f;
            }

            /**
             * @return First timer call since latest reset
             */
            double inline GetFirstCallTime()
            {
                return first_call_time_;
            }

            /**
             * @return Last timer call since latest reset
             */
            double inline GetLastCallTime()
            {
                return last_call_time_;
            }

            /**
             * @return Number of calls since lates reset
             */
            double inline GetRepetitionCount()
            {
                return repetition_count_;
            }

            /**
             * @return Sum of all calls since latest reset
             */
            double inline GetAggregatedTime()
            {
                return inner_calls_sum_;
            }

            /**
             * @return Duration between first and last call since latest reset
             */
            double inline GetTotalTime()
            {
                return last_call_time_ - first_call_time_;
            }

            /**
             * @return Average time over aggregated time
             */
            double inline GetAverageTimeOverAggregated()
            {
                return inner_calls_sum_/repetition_count_;
            }

            /**
             * @return Average time over total time
             */
            double inline GetAverageTimeOverTotal()
            {
                return GetTotalTime()/repetition_count_;
            }
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
             * Pompi timer
             * Used for all time measurements
             * @see Start()
             * @see End()
             */
            PompiTimer timer_;

            /**
             * Vector of papi event codes.
             * Dependant on PAPI_EVENTS enviroment variable.
             */
            std::vector< int > papi_events_;

            /**
             * Vector of papi event sets for each thread
             */
            std::vector< int > event_sets_;

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

            /**
             * Maximum number of events according to number of available hw counters
             */
            int max_events_;

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
             * Provides accumulated execution time betwwn all pairs of Start()
             * and Stop() methods. This behaviour is reset by ClearTimers().
             * @return execution time
             * @see Start()
             * @see Stop()
             * @see ClearTimers()
             */
            double inline GetExecutionTime();

            /**
             * Similar to GetExecutionTime. Used to calculate execution time 
             * averaged on number of trials.
             * @return        average execution time.
             * @see GetExecutionTime()
             */
            double inline GetAverageExecutionTime();

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
            void inline ClearTimers();

            /**
             * Prints results on standard output. Can be used to output counters
             * of any singular thread or all threads aggregated (implicit).
             * @param value         Any value selected by user (eg. total number of threads).
             * @param thread_id     ID of a thread to be outputted (invalid thread ID
             *                      results in aggregated output on all threads).
             */
            template < typename T >
            void PrintResults(T value, int thread_id = ALL_THREADS);

            /**
             * Simillar to PrintResults. Takes additional file description 
             * parameters to define output file.
             * @param value         Any value selected by user (eg. total number of threads).
             * @param file_name     C string name of a output file.
             * @param format        Output format.
             * @param thread_id     ID of a thread to be outputted (invalid thread ID
             *                      results in aggregated output on all threads).
             */
            template < typename T >
            void PrintResultsToFile(T value, const char * file_name, OutputFormat format, int thread_id = ALL_THREADS);

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
             * Private generic method that computes derived stat value.
             * @param  stat      Derived stat to be computed
             * @param  thread_id Derived stat will be computed for a thread with
             *                   this thread_id. In case thread_id is invalid,
             *                   the stat will be computed aggregated for all threads.
             * @return           Derived stat value
             */
            double ComputeDerivedStat(PapiDerivedStat stat, int thread_id = ALL_THREADS);

            /**
             * Private generic method defining print format for GNUPLOT
             * @param value         Any value selected by user (eg. total number of threads).
             * @param output        Output stream
             * @param thread_id     Output will be performed for thread with this ID.
             *                      In case thread_id is invalid, output will be
             *                      performed aggregated for all threads.
             */
            template < typename T >
            void PrintGnuplot(T value, std::ofstream &output, int thread_id = ALL_THREADS);
    };

    
    /////////////////////////////////////
    //  Implementation PUBLIC METHODS  //
    /////////////////////////////////////

    
    Base::Base()
    {
        max_threads_ = omp_get_max_threads();
        std::clog << "[Log] Maximum number of threads is " << max_threads_ << std::endl;

        thread_data_.resize(max_threads_);
        
        event_sets_.resize(max_threads_);
        for(int set = 0; set < event_sets_.size(); ++set)
            event_sets_[set] = PAPI_NULL;

        timer_.ResetTimer();

        int PAPI_return_value = PAPI_library_init(PAPI_VER_CURRENT);
        if(PAPI_return_value != PAPI_VER_CURRENT)
        {
            std::cerr << "[Error] Could not initialize PAPI library" << std::endl;
            exit(1);
        }
        std::clog << "[Log] PAPI library successfully initialized" << std::endl;

        max_events_ = PAPI_num_counters();
        std::clog << "[Log] Maximum number of counters is " << max_events_ << std::endl;


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
        if(papi_events_.size() > max_events_)
        {
            std::cerr << "[Warning] Cannot add any more events, skipping " << event << std::endl;
            return;
        }

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
        PAPI_register_thread();
        int thread_id = PAPI_thread_id();

        int PAPI_return_value = PAPI_create_eventset(&event_sets_[thread_id]);
        if(PAPI_return_value != PAPI_OK)
        {
            #pragma omp critical
            std::cerr << "[Error] Failed to initialize PAPI event set for thread #" << thread_id << std::endl;
        }

        PAPI_return_value = PAPI_add_events(event_sets_[thread_id], &papi_events_[0], papi_events_.size());
        if(PAPI_return_value != PAPI_OK)
        {
            #pragma omp critical
            std::cerr << "[Error] Failed to add events to event set for thread #" << thread_id << std::endl;
        }

        #pragma omp barrier

        #pragma omp single
        timer_.BeginTiming();

        PAPI_return_value = PAPI_start(event_sets_[thread_id]);
        if(PAPI_return_value != PAPI_OK)
        {
            #pragma omp critical
            std::cerr << "[Error] Failed to start counting in thread #" << thread_id << std::endl;
        }
    }


    void Base::Stop()
    {
        long long counter_values[papi_events_.size()];
        int thread_id = PAPI_thread_id();

        int PAPI_return_value = PAPI_stop(event_sets_[thread_id], counter_values);
        if(PAPI_return_value != PAPI_OK)
        {
            #pragma omp critical
            std::cerr << "[Error] Failed to stop counting in thread #" << thread_id << std::endl;
        }

        #pragma omp barrier

        #pragma omp single
        timer_.EndTiming();


        for(int event = 0; event < papi_events_.size(); ++event)
            thread_data_[thread_id][event] += counter_values[event];

        PAPI_return_value = PAPI_cleanup_eventset(event_sets_[thread_id]);
        if(PAPI_return_value != PAPI_OK)
        {
            #pragma omp critical
            std::cerr << "[Error] Failed to clean up event set in thread #" << thread_id << std::endl;
        }

        PAPI_return_value = PAPI_destroy_eventset(&event_sets_[thread_id]);
        if(PAPI_return_value != PAPI_OK)
        {
            #pragma omp critical
            std::cerr << "[Error] Failed to destroy event set in thread #" << thread_id << std::endl;
        }

        PAPI_unregister_thread();
    }


    double inline Base::GetExecutionTime()
    {
        return timer_.GetAggregatedTime();
    }


    double inline Base::GetAverageExecutionTime()
    {
        return timer_.GetAverageTimeOverAggregated();
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


    
    void inline Base::ClearTimers()
    {
        timer_.ResetTimer();
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
                return 100.0 * (counters[GetEventIndex(PAPI_L1_TCM)] / (double)(counters[GetEventIndex(PAPI_SR_INS)] + counters[GetEventIndex(PAPI_LD_INS)]));
            }
            case D_L2_TMR: {
                return 100.0 * (counters[GetEventIndex(PAPI_L2_TCM)] / (double)counters[GetEventIndex(PAPI_L2_TCA)]);
            }
            case D_L3_TMR: {
                return 100.0 * (counters[GetEventIndex(PAPI_L3_TCM)] / (double)counters[GetEventIndex(PAPI_L3_TCA)]);
            }
            case D_L1_DMR: {
                return 100.0 * (counters[GetEventIndex(PAPI_L1_DCM)] / (double)counters[GetEventIndex(PAPI_L1_DCA)]);
            }
            case D_L2_DMR: {
                return 100.0 * (counters[GetEventIndex(PAPI_L2_DCM)] / (double)counters[GetEventIndex(PAPI_L2_DCA)]);
            }
            case D_L3_DMR: {
                return 100.0 * (counters[GetEventIndex(PAPI_L3_DCM)] / (double)counters[GetEventIndex(PAPI_L3_DCA)]);
            }
            case D_BR_MPR: {
                if(EventAvailable(PAPI_BR_CN))
                    return 100.0 * (counters[GetEventIndex(PAPI_BR_MSP)] / (double)counters[GetEventIndex(PAPI_BR_CN)]);
                else
                    return 100.0 * (counters[GetEventIndex(PAPI_BR_MSP)] / (double)(counters[GetEventIndex(PAPI_BR_MSP)] + counters[GetEventIndex(PAPI_BR_PRC)]));
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

    template < typename T >
    void Base::PrintResults(T value, int thread_id)
    {
        std::vector< PapiDerivedStat > stats;
        GetDerivedStats(stats);

        if((thread_id >= 0)&&(thread_id <= max_threads_))
            std::cout << "Results for thread #" << thread_id << std::endl;
        else
        {
            std::cout << "Aggregated results on all threads" << std::endl;
            thread_id = ALL_THREADS;
        }

        long long results[papi_events_.size()];
        GetCounters(results, thread_id);

        std::cout << "Parameter value: " << value << std::endl;

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


    template < typename T >
    void Base::PrintResultsToFile(T value, const char * file_name, OutputFormat format, int thread_id)
    {
        if((thread_id < 0)||(thread_id > max_threads_))
            thread_id = ALL_THREADS;

        std::ofstream output;
        output.open(file_name, std::ofstream::app);

        // Add new output format methods here
        switch(format)
        {
            case GNUPLOT: {
                PrintGnuplot(value, output, thread_id);
                break;
            }
        }

        output.close();
    }


    // Format specific print methods
    template < typename T >
    void Base::PrintGnuplot(T value, std::ofstream &output, int thread_id)
    {
        std::vector< PapiDerivedStat > stats;
        GetDerivedStats(stats);

        output << '#';
        if(thread_id != ALL_THREADS)
            output << "THREAD";

        output << std::setw(15) << "VALUE";

        for(int event = 0; event < papi_events_.size(); ++event)
            output << std::setw(16) << papi_event_names_[event];
        for(int d_event = 0; d_event < stats.size(); ++d_event)
            output << std::setw(16) << GetDerivedStatName(stats[d_event]);
        output << std::endl;

        long long results[papi_events_.size()];
        GetCounters(results, thread_id);

        if(thread_id != ALL_THREADS)
            output << std::setw(7) << thread_id;

        output << std::setw(16) << value;

        for(int event = 0; event < papi_events_.size(); ++event)
            output << std::setw(16) << results[event];
        for(int d_event = 0; d_event < stats.size(); ++d_event)
            output << std::setw(16) << ComputeDerivedStat(stats[d_event], thread_id);
        output << std::endl;
    }
}
