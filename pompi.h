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
#include <vector>
#include <string>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#ifdef _OPENMP
  #include <omp.h>
#endif

/// pompi namespace.
namespace pompi
{

    /// Used in methods that support output on either only one thread or all threads.
    #define ALL_THREADS -1

    // Arbitrary boundary for max events
    #define MAX_EVENTS  20

    int g_eventset = PAPI_NULL;
#ifdef _OPENMP
    #pragma omp threadprivate(g_eventset)
#endif

    long long g_start_counter_values[MAX_EVENTS];
#ifdef _OPENMP
    #pragma omp threadprivate(g_start_counter_values)
#endif

    long long g_end_counter_values[MAX_EVENTS];
#ifdef _OPENMP
    #pragma omp threadprivate(g_start_counter_values)
#endif

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
        D_GIPC,     /**< Derived graduated instructions per second*/
        D_IIPC,     /**< Derived issued instructions per second*/
        D_BR_MPR,   /**< Derived branch missprediction rate [%]*/
        D_MFLOPS,   /**< Derived MFLOPS*/
        D_GFLOPS,   /**< Derived GFLOPS*/
        D_MIPS,     /**< Derived MIPS*/
        D_GIPS,     /**< Derived GIPS*/
        D_TLBM_PC,  /**< Derived TLB data misses per count */
        D_TLBM_PS,  /**< Derived TLB data misses per second */
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
            long long first_call_time_usec_;

            /**
             * Last call to timer
             */
            double last_call_time_;
            long long last_call_time_usec_;

            /**
             * Inner timer call (does not change first call)
             */
            double inner_call_time_;
            long long inner_call_time_usec_;

            /**
             * Inner calls aggregated duration
             */
            double inner_calls_sum_;
            long long inner_calls_sum_usec_;

            /**
             * Inner calls count
             */
            unsigned long int repetition_count_;


        public:

            /**
             * A constructor
             * Initializes all counters to zero
             */
            PompiTimer()
            {
                first_call_time_ = 0.0f;
                first_call_time_usec_ = 0;
                last_call_time_ = 0.0f;
                last_call_time_usec_ = 0;
                inner_call_time_ = 0.0f;
                inner_call_time_usec_ = 0;
                inner_calls_sum_ = 0.0f;
                inner_calls_sum_usec_ = 0;
                repetition_count_ = 0;
            }

            /**
             * Begins timing
             * Sets first call if it is zero, otherwise sets inner call time
             */
            void BeginTiming()
            {
#ifdef _OPENMP
                inner_call_time_ = omp_get_wtime();
                if(first_call_time_ == 0.0f)
                    first_call_time_ = inner_call_time_;
#endif
                inner_call_time_usec_ = PAPI_get_real_usec();
                if(first_call_time_usec_ == 0)
                    first_call_time_usec_ = inner_call_time_usec_;
            }

            /**
             * Ends timing
             */
            void EndTiming()
            {
#ifdef _OPENMP
                last_call_time_ = omp_get_wtime();
                inner_calls_sum_ += last_call_time_ - inner_call_time_;
#endif
                last_call_time_usec_ = PAPI_get_real_usec();
                inner_calls_sum_usec_ += last_call_time_usec_ - inner_call_time_usec_;
                repetition_count_++;
            }

            /**
             * Clears entire timer
             */
            void ResetTimer()
            {
                first_call_time_ = 0.0f;
                first_call_time_usec_ = 0;
                last_call_time_ = 0.0f;
                last_call_time_usec_ = 0;
                inner_call_time_ = 0.0f;
                inner_call_time_usec_ = 0;
                inner_calls_sum_ = 0.0f;
                inner_calls_sum_usec_ = 0;
                repetition_count_ = 0;
            }

            /**
             * @return First timer call since latest reset
             */
            double GetFirstCallTime()
            {
#ifndef _OPENMP
                return first_call_time_usec_ / 1e6;
#endif
                return first_call_time_;
            }

            long long GetFirstCallTimeUsec()
            {
                return first_call_time_usec_;
            }

            /**
             * @return Last timer call since latest reset
             */
            double GetLastCallTime()
            {
#ifndef _OPENMP
                return last_call_time_usec_ / 1e6;
#endif
                return last_call_time_;
            }

            long long GetLastCallTimeUsec()
            {
                return last_call_time_usec_;
            }

            /**
             * @return Number of calls since lates reset
             */
            double GetRepetitionCount()
            {
                return repetition_count_;
            }

            void SetRepetitionCount(unsigned long int count)
            {
                repetition_count_ = count;
            }

            /**
             * @return Sum of all calls since latest reset
             */
            double GetAggregatedTime()
            {
#ifndef _OPENMP
                return inner_calls_sum_usec_ / 1e6;
#endif
                return inner_calls_sum_;
            }

            long long GetAggregatedTimeUsec()
            {
                return inner_calls_sum_usec_;
            }

            /**
             * @return Duration between first and last call since latest reset
             */
            double GetTotalTime()
            {
#ifndef _OPENMP
                return GetTotalTimeUsec() / 1e6;
#endif
                return last_call_time_ - first_call_time_;
            }

            long long GetTotalTimeUsec()
            {
                return last_call_time_usec_ - first_call_time_usec_;
            }

            /**
             * @return Average time over aggregated time
             */
            double GetAverageTimeOverAggregated()
            {
#ifndef _OPENMP
                return GetAverageTimeOverAggregatedUsec() / 1e6;
#endif
                return inner_calls_sum_/repetition_count_;
            }

            double GetAverageTimeOverAggregatedUsec()
            {
                return inner_calls_sum_usec_/repetition_count_;
            }

            /**
             * @return Average time over total time
             */
            double GetAverageTimeOverTotal()
            {
#ifndef _OPENMP
                return GetAverageTimeOverTotalUsec() / 1e6;
#endif
                return GetTotalTime()/repetition_count_;
            }

            double GetAverageTimeOverTotalUsec()
            {
                return GetTotalTimeUsec()/repetition_count_;
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
            int hw_counters_;


        public:

            /**
             * Base constructor
             * Initializes PAPI library and PAPI library threading
             * support for OpenMP. Parses PAPI_EVENTS enviroment variable
             * and modifies papi_events_ and papi_event_names_ accordingly.
             */
            Base();

            /**
             * Base destructor
             */
            ~Base()
            {
                PAPI_shutdown();
            }

        public:

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
             * Sets repetition count
             * @param count new count
             */
            void SetRepetitionCount(unsigned long int count);

            /**
             * Provides accumulated execution time betwwn all pairs of Start()
             * and Stop() methods. This behaviour is reset by ClearTimers().
             * @return execution time
             * @see Start()
             * @see Stop()
             * @see ClearTimers()
             */
            double GetExecutionTime();

            /**
             * Similar to GetExecutionTime. Used to calculate execution time
             * averaged on number of trials.
             * @return        average execution time.
             * @see GetExecutionTime()
             */
            double GetAverageExecutionTime();

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
            void PrintGnuplot(T value, FILE * output, int thread_id = ALL_THREADS);
    };


    /////////////////////////////////////
    //  Implementation PUBLIC METHODS  //
    /////////////////////////////////////


    Base::Base()
    {

#ifdef _OPENMP
        max_threads_ = 24;//omp_get_max_threads();
        printf("[Log] OpenMP has been found, Maximum number of threads is %d\n", max_threads_);
#else
        max_threads_ = 1;
        printf("[Log] OpenMP not found, POMPI will be operating in single thread mode\n");
#endif
        thread_data_.resize(max_threads_);

        timer_.ResetTimer();

        int PAPI_return_value = PAPI_library_init(PAPI_VER_CURRENT);
        if(PAPI_return_value != PAPI_VER_CURRENT)
        {
            fprintf(stderr, "[Error] Could not initialize PAPI library\n");
            exit(1);
        }
        printf("[Log] PAPI library successfully initialized\n");

        hw_counters_ = PAPI_num_counters();
        if(hw_counters_ > MAX_EVENTS)
            hw_counters_ = MAX_EVENTS;

        printf("[Log] Maximum number of events is %d\n", hw_counters_);

#ifdef _OPENMP
        PAPI_return_value = PAPI_thread_init((unsigned long (*)(void)) (omp_get_thread_num));
        if(PAPI_return_value != PAPI_OK)
        {
            fprintf(stderr, "[Error] Coult not initialize OpenMP threading for PAPI\n");
            exit(1);
        }
        printf("[Log] OpenMP threading for PAPI successfully initialized\n");
#endif
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
        if(papi_events_.size() > hw_counters_)
        {
            fprintf(stderr, "[Warning] Cannot add any more events, skipping %s\n", event);
            return;
        }

        int event_id;
        int PAPI_return_value = PAPI_event_name_to_code(event, &event_id);
        if(PAPI_return_value != PAPI_OK)
            fprintf(stderr, "[Warning] Papi event `%s` does not exists, skipping\n", event);
        else
        {
            if(std::find(papi_events_.begin(), papi_events_.end(), event_id) == papi_events_.end())
            {
                papi_event_names_.push_back(std::string(event));
                papi_events_.push_back(event_id);
                for(int i = 0; i < max_threads_; ++i)
                    thread_data_[i].push_back(0);
                printf("[Log] Adding papi event `%s`\n", event);
            }
            else
                fprintf(stderr, "[Warning] Papi event `%s` already listed, skipping\n", event);
        }
    }


    void Base::Start()
    {
        timer_.BeginTiming();

#ifdef _OPENMP
        #pragma omp parallel
        {
            PAPI_register_thread();
            int thread_id = PAPI_thread_id();

            int PAPI_return_value = PAPI_create_eventset(&g_eventset);
            if(PAPI_return_value != PAPI_OK)
            {
                #pragma omp critical
                fprintf(stderr, "[Error] Failed to initialize PAPI event set for thread #%d\n", thread_id);
            }

            PAPI_return_value = PAPI_add_events(g_eventset, &papi_events_[0], papi_events_.size());
            if(PAPI_return_value != PAPI_OK)
            {
                #pragma omp critical
                fprintf(stderr, "[Error] Failed to add events to event set for thread #%d\n", thread_id);
            }

            PAPI_return_value = PAPI_start(g_eventset);
            if(PAPI_return_value != PAPI_OK)
            {
                #pragma omp critical
                fprintf(stderr, "[Error] Failed to start counting in thread #%d\n", thread_id);
            }

            PAPI_return_value = PAPI_read(g_eventset, g_start_counter_values);
            if(PAPI_return_value != PAPI_OK)
            {
                #pragma omp critical
                fprintf(stderr, "[Error] Failed to read counters in thread #%d\n", thread_id);
            }
        }
#else
        int PAPI_return_value = PAPI_create_eventset(&g_eventset);
        if(PAPI_return_value != PAPI_OK)
            fprintf(stderr, "[Error] Failed to initialize PAPI event set\n");

        PAPI_return_value = PAPI_add_events(g_eventset, &papi_events_[0], papi_events_.size());
        if(PAPI_return_value != PAPI_OK)
            fprintf(stderr, "[Error] Failed to add events to event set\n");

        PAPI_return_value = PAPI_start(g_eventset);
        if(PAPI_return_value != PAPI_OK)
            fprintf(stderr, "[Error] Failed to start counting\n");

        PAPI_return_value = PAPI_read(g_eventset, g_start_counter_values);
        if(PAPI_return_value != PAPI_OK)
            fprintf(stderr, "[Error] Failed to read counters\n");
#endif
    }


    void Base::Stop()
    {
#ifdef _OPENMP
        #pragma omp parallel
        {
            int thread_id = PAPI_thread_id();

            int PAPI_return_value = PAPI_read(g_eventset, g_end_counter_values);
            if(PAPI_return_value != PAPI_OK)
            {
                #pragma omp critical
                fprintf(stderr, "[Error] Failed to read counters in thread #%d\n", thread_id);
            }

            for(int event = 0; event < papi_events_.size(); ++event)
                thread_data_[thread_id][event] += g_end_counter_values[event];

            PAPI_return_value = PAPI_stop(g_eventset, NULL);
            if(PAPI_return_value != PAPI_OK)
            {
                #pragma omp critical
                fprintf(stderr, "[Error] Failed to stop counting in thread #%d\n", thread_id);
            }

            PAPI_return_value = PAPI_cleanup_eventset(g_eventset);
            if(PAPI_return_value != PAPI_OK)
            {
                #pragma omp critical
                fprintf(stderr, "[Error] Failed to clean up event set in thread #%d\n", thread_id);
            }

            PAPI_return_value = PAPI_destroy_eventset(&g_eventset);
            if(PAPI_return_value != PAPI_OK)
            {
                #pragma omp critical
                fprintf(stderr, "[Error] Failed to destroy event set in thread #%d\n", thread_id);
            }

            PAPI_unregister_thread();
        }
#else
        int PAPI_return_value = PAPI_read(g_eventset, g_end_counter_values);
        if(PAPI_return_value != PAPI_OK)
            fprintf(stderr, "[Error] Failed to read counters\n");

        for(int event = 0; event < papi_events_.size(); ++event)
            thread_data_[0][event] += g_end_counter_values[event];

        PAPI_return_value = PAPI_stop(g_eventset, NULL);
        if(PAPI_return_value != PAPI_OK)
            fprintf(stderr, "[Error] Failed to stop counting\n");

        PAPI_return_value = PAPI_cleanup_eventset(g_eventset);
        if(PAPI_return_value != PAPI_OK)
            fprintf(stderr, "[Error] Failed to clean up event set\n");

        PAPI_return_value = PAPI_destroy_eventset(&g_eventset);
        if(PAPI_return_value != PAPI_OK)
            fprintf(stderr, "[Error] Failed to destroy event set\n");
#endif
        timer_.EndTiming();
    }


    void Base::SetRepetitionCount(unsigned long int count)
    {
        timer_.SetRepetitionCount(count);
    }


    double Base::GetExecutionTime()
    {
        return timer_.GetAggregatedTime();
    }


    double Base::GetAverageExecutionTime()
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


    void  Base::ClearTimers()
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

        if (EventAvailable(PAPI_TOT_INS) && EventAvailable(PAPI_TOT_CYC))
            stats.push_back(D_GIPC);

        if (EventAvailable(PAPI_TOT_IIS) && EventAvailable(PAPI_TOT_CYC))
            stats.push_back(D_IIPC);

        if (EventAvailable(PAPI_FP_INS))
        {
            stats.push_back(D_MFLOPS);
            stats.push_back(D_GFLOPS);
        }

        if (EventAvailable(PAPI_TOT_INS))
        {
            stats.push_back(D_MIPS);
            stats.push_back(D_GIPS);
        }

        if (EventAvailable(PAPI_TLB_DM))
        {
            stats.push_back(D_TLBM_PC);
            stats.push_back(D_TLBM_PS);
        }
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
            case D_GIPC: {
                return (static_cast<double>(counters[GetEventIndex(PAPI_TOT_INS)]) / counters[GetEventIndex(PAPI_TOT_CYC)]);
            }
            case D_IIPC: {
                return (static_cast<double>(counters[GetEventIndex(PAPI_TOT_IIS)]) / counters[GetEventIndex(PAPI_TOT_CYC)]);
            }
            case D_MFLOPS: {
                return (static_cast<double>(counters[GetEventIndex(PAPI_FP_INS)]) / timer_.GetAggregatedTime()) / 1e6;
            }
            case D_GFLOPS: {
                return (static_cast<double>(counters[GetEventIndex(PAPI_FP_INS)]) / timer_.GetAggregatedTime()) / 1e9;
            }
            case D_MIPS: {
                return (static_cast<double>(counters[GetEventIndex(PAPI_TOT_INS)]) / timer_.GetAggregatedTime()) / 1e6;
            }
            case D_GIPS: {
                return (static_cast<double>(counters[GetEventIndex(PAPI_TOT_INS)]) / timer_.GetAggregatedTime()) / 1e9;
            }
            case D_TLBM_PC: {
                return (static_cast<double>(counters[GetEventIndex(PAPI_TLB_DM)]) / timer_.GetRepetitionCount());
            }
            case D_TLBM_PS: {
                return (static_cast<double>(counters[GetEventIndex(PAPI_TLB_DM)]) / timer_.GetAggregatedTime());
            }
        }

        return 0;
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
            case D_GIPC: {
                return "D_GIPC";
            }
            case D_IIPC: {
                return "D_IIPC";
            }
            case D_MFLOPS: {
                return "D_MFLOPS";
            }
            case D_GFLOPS: {
                return "D_GFLOPS";
            }
            case D_MIPS: {
                return "D_MIPS";
            }
            case D_GIPS: {
                return "D_GIPS";
            }
            case D_TLBM_PC: {
                return "D_TLBM_PC";
            }
            case D_TLBM_PS: {
                return "D_TLBM_PS";
            }
        }
        return "";
    }

    template < typename T >
    void Base::PrintResults(T value, int thread_id)
    {
        std::vector< PapiDerivedStat > stats;
        GetDerivedStats(stats);

        if((thread_id >= 0)&&(thread_id <= max_threads_))
            printf("Results for thread #%d\n", thread_id);
        else
        {
            printf("Aggregated results on all threads\n");
            thread_id = ALL_THREADS;
        }

        double t_value = static_cast<double>(value);
        printf("Parameter value: %f\n", t_value);

        long long results[papi_events_.size()];
        GetCounters(results, thread_id);

        for(int event = 0; event < papi_events_.size(); ++event)
            printf("%-15s:%20d\n", papi_event_names_[event].c_str(), results[event]);

        for(int d_event = 0; d_event < stats.size(); ++d_event)
            printf("%-15s:%20f\n", GetDerivedStatName(stats[d_event]).c_str(), ComputeDerivedStat(stats[d_event], thread_id));
    }


    template < typename T >
    void Base::PrintResultsToFile(T value, const char * file_name, OutputFormat format, int thread_id)
    {
        if((thread_id < 0)||(thread_id > max_threads_))
            thread_id = ALL_THREADS;

        FILE * output;
        output = fopen(file_name, "a");

        // Add new output format methods here
        switch(format)
        {
            case GNUPLOT: {
                PrintGnuplot(value, output, thread_id);
                break;
            }
        }

        fclose(output);
    }


    // Format specific print methods
    template < typename T >
    void Base::PrintGnuplot(T value, FILE * output, int thread_id)
    {
        std::vector< PapiDerivedStat > stats;
        GetDerivedStats(stats);

        fprintf(output, "#");
        if(thread_id != ALL_THREADS)
            fprintf(output, "THREAD");

        fprintf(output, "%15s", "VALUE");

        for(int event = 0; event < papi_events_.size(); ++event)
            fprintf(output, "%20s", papi_event_names_[event].c_str());
        for(int d_event = 0; d_event < stats.size(); ++d_event)
            fprintf(output, "%20s", GetDerivedStatName(stats[d_event]).c_str());
        fprintf(output, "\n");

        long long results[papi_events_.size()];
        GetCounters(results, thread_id);

        if(thread_id != ALL_THREADS)
            fprintf(output, "%-7d", thread_id);

        double t_value = static_cast<double>(value);
        fprintf(output, "%16f", t_value);

        for(int event = 0; event < papi_events_.size(); ++event)
            fprintf(output, "%20d", results[event]);
        for(int d_event = 0; d_event < stats.size(); ++d_event)
            fprintf(output, "%20f", ComputeDerivedStat(stats[d_event], thread_id));
        fprintf(output, "\n");
    }
}
