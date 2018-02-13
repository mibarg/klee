//
// Created by liav on 12/02/18.
//

#ifndef KLEE_STRINGUTILS_H
#define KLEE_STRINGUTILS_H

/* String-related methods */

/**
Count occurences of substr in str
*/
static int count_occurences(const char * str, char * substr);

/**
Count the number of capital letters in a string
*/
static int count_capitals(const char * str) ;

/**
Given a String str, a search_char and a break char, calculates the
length of all possible sequences of consequtive search_char's without a middle break_char.
These possible sequences of length are aggregated into the std, avg and max length,
and saved into the first three places in results

@assert len(results) >= 3
*/
static void consecutive_occurences_stats(const char * str, char search_char, char break_char, double * results);
#endif //KLEE_STRINGUTILS_H
