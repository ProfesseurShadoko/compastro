#pragma once

#include <string> 
#include <iostream>
#include <chrono>
#include <thread>

/**
 * @brief Class to color a string in the terminal. Call ColoredString("hop")::green()
 * and then print it to get a green "hop" in the terminal. All functions return type string.
 * 
 * 
*/
class ColoredString {
    private:
        std::string str;
        static const std::string _green;
        static const std::string _red;
        static const std::string _yellow;
        static const std::string _blue;
        static const std::string _purple;
        static const std::string _cyan;
        static const std::string _reset;

    public:
        ColoredString(const std::string str);

        std::string green();
        std::string red();
        std::string yellow();
        std::string blue();
        std::string purple();
        std::string cyan();
};


ColoredString cstr(std::string str);
ColoredString cstr(int i);
ColoredString cstr(double d);
ColoredString cstr(ColoredString cs);


class MutableClass {
    private:
        static int mute_count;

    public:
        /**
         * @brief Increment the mute_count static variable. This prevents the message from being printed.
        */
        static void mute();

        /**
         * @brief Decrement the mute_count static variable. This allows the message to be printed,
         * provided mute() method has not been called more times than unmute() method.
        */
        static void unmute();

        /**
         * @brief Check if the message is muted by looking at the private static variable mute_count.
         * This variable gets incremented or decremented by the mute() and unmute() methods.
         * @return true if the message is muted, false otherwise
        */
        static bool is_muted();

        /**
         * @brief Create a ColoredString object from a string.
         * @return ColoredString object. This means you can use:
         * MutableCLass::print(MutableClass::cstr("hello")->green())
        */
        static ColoredString cstr(std::string str);

        /**
         * @brief Print the message if it is not muted.
        */
        static void print(std::string msg = "");
};


class Message : public MutableClass {
private:
    std::string msg;
    std::string typ;

public:
    Message(const std::string msg, std::string typ = "i");
};



class Timer {

private:
    std::string name;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_time;
    std::string status;

public:
    Timer();
    Timer(std::string name);
    void start();
    void stop();
    void display();

    /**
     * @brief Get the time in nanoseconds between the start and stop of the timer.
     */
    int getTime();
};



class ProgressBar {
private:
    int length;
    int progress;
    Timer timer;
    void display();
    static const int bar_length = 50;
    static bool muted;


public:
    ProgressBar(int length);
    void update();
    static void sleep(int ms);
    static void mute();
    static void unmute();
};