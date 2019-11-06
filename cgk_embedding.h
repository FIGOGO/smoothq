#pragma once
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <algorithm>
#include "util.h"


std::vector<int> generate_random_sampling_positions(int size_cgk, int size_smq)
{
    std::vector<int> hash_lsh;
    hash_lsh.reserve(size_cgk);
    for (int j = 0; j < size_cgk; j++)
        hash_lsh.push_back(j);
    random_shuffle(hash_lsh.begin(), hash_lsh.end());
    hash_lsh.resize(size_smq);
    sort(hash_lsh.begin(), hash_lsh.end());

    return hash_lsh;
};

int** generate_random_binary(int size_cgk) {
    const int num_char = 5;
    int **p = new int *[num_char];
    for (int t = 0; t < num_char; t++) {
        p[t] = new int [size_cgk];
        for (int d = 0; d < size_cgk; d++) {
            p[t][d] = rand()%2;
        }
    }
    return p;
};

std::string cgk_embedding (const std::string& input, int size_cgk, int** random_string) {
    std::string output;
    for (int i = 0, j = 0; j < size_cgk; j++) {
        char s = i < input.size() ? input.at(i) : 'N';
        i = i + *(random_string[dict[s]] + j);
        output += s;
    }
    return output;
};

std::string cgk_sampling (const std::string& input, std::vector<int>& random_positions) {
    std::string output;
    for (int index: random_positions) output += input.at(index);
    return output;
};


