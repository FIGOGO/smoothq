#pragma once

#include <string.h>
#include <algorithm>
#include <omp.h>

typedef int32_t int32;

// For two strings aligned from the beginning, find the index of their first difference
int find_first_diff(const char *x, const int x_len, const  char *y, const int y_len) {
    // the procedure assumes little-endian memory layout for integers (intel x86 type).
    // For big-endian one would need to modify the computation using __builtin_ctz to something using __builtin_clz
    int min_len = std::min(x_len, y_len);
    int steps = min_len / 4;
    int remainder = min_len % 4;

    int index = 0;
    for (int i = 0; i < steps; i++) {
        if (*((int32*)x) == *((int32*)y)) {
            x += 4; y += 4; index += 4;    // find the first 8-bytes that differ
        } else {
            break;
        }
    }

    if (index / 4 < steps) {  // Not all first 4*steps letters are the same
        index += (__builtin_ctz(*((int32*)x) - *((int32*)y)) >> 3);
        return index;
    }

    for (int i = 0; i < remainder; i++) {  // Continue on the remainder letters (at most 3 letters)
        if (*x == *y) {
            x++;  y++;  index++;  // only move by 1
        } else {
            break;
        }
    }

    return index;
}

bool edit_le_eq_one(const char *x, const int x_len, const  char *y, const int y_len) {
    // Make y_len the larger one
    if (x_len > y_len) return edit_le_eq_one(y, y_len, x, x_len);

    int len_diff = y_len - x_len;
    if (len_diff > 1) return false;

    int first_diff = find_first_diff(x, x_len, y, y_len);

    // Equal length case
    if (len_diff == 0) {
        if (first_diff >= x_len-1) return true;
        int second_diff = find_first_diff(x+first_diff+1, x_len-first_diff-1, y+first_diff+1, y_len-first_diff-1);
        int third_diff = find_first_diff(x+first_diff+second_diff+1, x_len-first_diff-second_diff-1,
                                         y+first_diff+second_diff+1, y_len-first_diff-second_diff-1);
        return (third_diff == x_len-first_diff-second_diff-1);
    }

    // Length difference = 1 case
    if (len_diff == 1) {
        if (first_diff == y_len) return true;
        int second_diff = find_first_diff(x+first_diff, x_len-first_diff, y+first_diff+1, y_len-first_diff-1);
        return (second_diff == x_len-first_diff);
    }

    return false;
}

// Check if two equal length string has edit distance less or equal 2
//  if true, will return 0
//  if false, will return -1
int edit_early_quit(const char *x, const char *y, const int len, int k=2) {
    int first_diff = find_first_diff(x, len, y, len);
    if (first_diff >= len-2) return 0;

    int new_len = len - first_diff;
    if (edit_le_eq_one(x+1+first_diff, new_len-1, y+1+first_diff, new_len-1) ||
        edit_le_eq_one(x+first_diff, new_len, y+1+first_diff, new_len-1) ||
        edit_le_eq_one(x+1+first_diff, new_len-1, y+first_diff, new_len))
        return 0;

    return -1;
}
