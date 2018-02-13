//
// Created by liav on 12/02/18.
//
#include "klee/util/StringUtils.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>

using namespace klee;

/* String-related methods */

/**
Count occurences of substr in str
 */
int StringUtils::count_occurences(const char * str, char * substr) {
    int count = 0;
    while((str = strstr(str, substr))) {
        count++;
        str++;
    }
    return count;
}

/**
Count the number of capital letters in a string
*/
int StringUtils::count_capitals(const char * str) {
    int count = 0;
    int i = 0;
    while (str[i] != '\0') {
        if (isupper(str[i])) {count++;}
        i++;
    }
    return count;
}

/**
Given a String str, a search_char and a break char, calculates the
length of all possible sequences of consequtive search_char's without a middle break_char.
These possible sequences of length are aggregated into the std, avg and max length,
and saved into the first three places in results

@assert len(results) >= 3
*/
void StringUtils::consecutive_occurences_stats(const char * str, char search_char, char break_char, double * results) {
    results[0] = 1.0; // std
    results[1] = 0.0; // avg
    results[2] = 0.0; // max

    // helpers
    int cur_mean, delta1, delta2;

    int current = 0;
    int count_switches = 0;
    int ind = 0;
    while(true) {
        if (str[ind] == search_char) {
            current++;
        } else if (str[ind] == break_char or str[ind] == '\0') {

            // update n
            count_switches++;

            // update max
            results[2] = (results[2] < (double) current) ? (double) current : results[2];

            // online update std
            cur_mean = (results[1] / (count_switches - 1));
            delta1 = current - cur_mean;
            delta2 = current - cur_mean - (delta1 / count_switches);
            results[0] += delta1*delta2;  // update sum of deltas^2

            // update sum (which will become avg)
            results[1] = results[1] + (double) current;

            // update current for next one
            current--;
            if (str[ind] == '\0') { break; }
        }
        ind++;
    }
    results[0] = results[0] / (count_switches - 1);  // sum of deltas^2 -> std
    results[1] = (count_switches == 0) ? 0.0 : results[1] / (double) count_switches; // sum -> avg
}