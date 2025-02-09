#include "message.hpp" // everything is declared there
#include <string> 
#include <iostream>
#include <sys/resource.h>

/**
 * Colored String
*/

const std::string ColoredString::_red = "\033[91m";
const std::string ColoredString::_green = "\033[92m";
const std::string ColoredString::_yellow = "\033[93m";
const std::string ColoredString::_blue = "\033[94m";
const std::string ColoredString::_purple = "\033[95m";
const std::string ColoredString::_cyan = "\033[96m";
const std::string ColoredString::_reset = "\033[0m";

ColoredString::ColoredString(const std::string str) {
    this->str = str;
}

std::string ColoredString::green() {
    return _green + str + _reset;
}

std::string ColoredString::red() {
    return _red + str + _reset;
}

std::string ColoredString::yellow() {
    return _yellow + str + _reset;
}

std::string ColoredString::blue() {
    return _blue + str + _reset;
}

std::string ColoredString::purple() {
    return _purple + str + _reset;
}

std::string ColoredString::cyan() {
    return _cyan + str + _reset;
}


ColoredString cstr(std::string str) {
    return ColoredString(str);
}

ColoredString cstr(int i) {
    return ColoredString(std::to_string(i));
}

ColoredString cstr(double d) {
    return ColoredString(std::to_string(d));
}

ColoredString cstr(ColoredString cs) {
    return cs;
}




/**
 * Mutable Class
*/

int MutableClass::mute_count = 0;

void MutableClass::mute() {
    mute_count++;
}

void MutableClass::unmute() {
    mute_count--;
    mute_count = mute_count < 0 ? 0 : mute_count;
}

bool MutableClass::is_muted() {
    return mute_count > 0;
}

ColoredString MutableClass::cstr(std::string str) {
    return ColoredString(str);
}

void MutableClass::print(std::string msg) {
    if (!is_muted()) {
        std::cout << msg << std::endl;
    }
}

/**
 * Message!
*/

Message::Message(const std::string msg, std::string typ) {
    this->msg = msg;
    std::string prefix;

    if (typ == "i") {
        prefix = cstr("[i] ").cyan();
    } else if (typ == "#") {
        prefix = cstr("[#] ").green();
    } else if (typ == "!") {
        prefix = cstr("[!] ").red();
    } else if (typ == "?") {
        prefix = cstr("[?] ").yellow();
    } else {
        throw std::invalid_argument("Invalid typ argument. Must be Must be one of '#', '!', '?', 'i'.");
    }

    print(prefix + this->msg);
}


/*
#############
### TIMER ###
#############
*/

Timer::Timer() : name("DEFAULT"), status("EMPTY") {};

Timer::Timer(std::string name) : name(name), status("EMPTY") {};

void Timer::start() {
    if (status == "RUNNING") {
        throw std::runtime_error("Timer is already running!");
    }
    status = "RUNNING";
    start_time = std::chrono::high_resolution_clock::now(); // reset here!
}

long long Timer::stop() {
    end_time = std::chrono::high_resolution_clock::now();
    if (status != "RUNNING") {
        throw std::runtime_error("Timer has not been started yet.");
    }
    status = "STOPPED";
    return getTimeNs();
}

void Timer::display() {
    if (status != "STOPPED") {
        throw std::runtime_error("Timer has not been stopped yet.");
    }
    long long ns_time = getTimeNs();
    long long ms_time = ns_time / 1e6; // static cast to int
    long long s_time = ns_time / 1e9;

    std::string display_name = (name=="DEFAULT" ? "" : "[" + name + "]");
    
    std::string prefix = Message::cstr("Timer ").purple() + display_name + ": ";

    if (s_time > 1) {
        Message::print(prefix + std::to_string(s_time) + "s");
    } else if (ms_time > 1) {
        Message::print(prefix + std::to_string(ms_time) + "ms");
    } else {
        Message::print(prefix + std::to_string(ns_time) + "ns");
    }
}

long long Timer::getTimeNs() {
    if (status != "STOPPED") {
        throw std::runtime_error("Timer has not been stopped yet.");
    }
    // time in ns
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
}

/**
 * ###################
 * ### PROGRESSBAR ###
 * ###################
 */

bool ProgressBar::muted = false;

ProgressBar::ProgressBar(int length) : length(length), progress(0) {
    display();
}

void ProgressBar::display() {

    if (muted) {
        return; // might be redundant with other stuff, but just in case
    }
    
    std::string incomplete = "━";
    std::string complete = cstr("━").blue();
    int n_complete = ((progress * bar_length) / length);

    std::string bar = "";
    // fill the bar with complete times █
    for (int i = 0; i < n_complete; i++) {
        bar += complete;
    }
    // fill the bar with incomplete times █
    for (int i = n_complete; i < bar_length; i++) {
        bar += incomplete;
    }

    int progress_percent = ((progress * 100) / length);
    if (progress >= length) {
        progress_percent = 100;
    }

    std::string progress_percent_str = cstr(std::to_string(progress_percent) + "%").red();
    std::string next_print = "\r" + cstr("[%]").blue() + " Progress: " + bar + " (" + progress_percent_str + ")";

    if (next_print != previous_print) {
        previous_print = next_print;

        // let's add some information about the memory usage
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        std::string memory_usage = " > " + std::to_string(usage.ru_maxrss / 1000) + " MB";

        std::cout << next_print + memory_usage << std::flush;

    }


    
}

void ProgressBar::update() {
    progress++;

    if (muted) {
        return; // no display, but still update
    }

    display();

    if (progress == length) {
        std::cout << std::endl;
    }
}



void ProgressBar::sleep(int ms) {
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

void ProgressBar::mute() {
    muted = true;
}

void ProgressBar::unmute() {
    muted = false;
}

void ProgressBar::print(std::string msg) {
   
    if (muted) {
        std::cout << msg << std::endl;
        return;
    }

    // print one empty line to erase the line. then go back to front of line.
    std::string empty_line = "\r";
    for (size_t i=0; i<bar_length*2; i++) {
        empty_line += " ";
    }
    std::cout << empty_line << std::flush;
    std::cout << "\r" << msg << std::endl;
    std::cout << previous_print << std::flush;
    
}



void ProgressBar::printMemoryUsage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::string to_print = "Memory Usage: " + std::to_string(usage.ru_maxrss) + " KB";
    print(to_print);
}