#ifndef __timing_h__
#define __timing_h__

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

namespace timing {

    class Timer {
    public:
        // Constructor: Starts the timer
        Timer() : start_time(std::chrono::high_resolution_clock::now()) {}

        // Restart the timer
        void restart() {
            start_time = std::chrono::high_resolution_clock::now();
        }

        // Print the elapsed time with dynamic formatting
        void printElapsed(const std::string& label = "Time") const {
            auto stop_time = std::chrono::high_resolution_clock::now();
            auto duration = stop_time - start_time;

            // Convert to different time units
            auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
            auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
            auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration).count();
            auto hours = std::chrono::duration_cast<std::chrono::hours>(duration).count();

            // Print the elapsed time in the most appropriate format
            std::cout << "\n" << label << ": ";
            if (microseconds < 1000) {
                std::cout << microseconds << " microseconds." << std::endl;
            } else if (milliseconds < 1000) {
                std::cout << milliseconds << " milliseconds." << std::endl;
            } else if (seconds < 60) {
                std::cout << seconds << " seconds." << std::endl;
            } else if (minutes < 60) {
                std::cout << minutes << " minutes and " << (seconds % 60) << " seconds." << std::endl;
            } else {
                std::cout << hours << " hours, " << (minutes % 60) << " minutes, and " << (seconds % 60) << " seconds." << std::endl;
            }
        }

    private:
        std::chrono::high_resolution_clock::time_point start_time;
    };

} // namespace timing

#endif // __timing_h__