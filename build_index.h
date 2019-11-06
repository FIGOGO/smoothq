#pragma once

#include <unordered_set>
#include <fstream>
#include <cstring>
#include "util.h"
#include "cgk_embedding.h"
#include "xxhash.h"
#include <algorithm>

template <class InIt, class OutIt, class size_type>
void window(InIt b, InIt e, OutIt d, size_type len) {
    for (auto s = std::next(b, len); s != e; ++s, ++b, ++d)
        *d = { b, s };
}

std::vector<std::string> generate_qgram(std::string const &in, size_t q) {
    std::string input = in;
    std::vector<std::string> words;
    window(input.begin(), input.end(), back_inserter(words), q);
    return words;
}

inline std::string embedstr(std::string& input, int **p, std::vector<int>& hash_lsh)
{
    std::string output;
    int i = 0, partdigit = 0;

    for (int j = 0; partdigit < hash_lsh.size(); ++j)
    {
        char s = i < input.size() ? input[i] : 'N';
        //only record the bits used in LSH
        while (j == hash_lsh[partdigit])
        {
            output += s;
            partdigit += 1;
        }
        i = i + *(p[dict[s]] + j);
    }
    return output;
}

std::vector<Signature_inter> generate_ith_signature(int i, int **p, std::vector<int>& hash_lsh, int seed) {
    std::vector<std::string> qgrams = generate_qgram(oridata[i].str, SIZE_Q);
    std::vector<Signature_inter> sig3;
    sig3.reserve(qgrams.size());
    for (int j = 0; j < qgrams.size(); ++j)
    {
        std::string smooth_qgram = embedstr(qgrams[j], p, hash_lsh);
        HASH_TYPE hash;
        const char *cstr = smooth_qgram.c_str();
        hash = XXH64(cstr, SIZE_SMQ, seed);
        sig3.push_back({i, j, hash});
    }
    qgrams.clear();

    std::sort(sig3.begin(), sig3.end(), sortSiginter);
    // Only keep a portion of signatures
    // Include the last hash values that are just beyond the threshold
    int num_hash = (int) ((oridata[i].len - SIZE_Q + 1) * SAMPLING_RATE) - 1;
    HASH_TYPE current_hash = sig3[num_hash].hash_value;
    while (num_hash < qgrams.size()
           && current_hash == sig3[num_hash].hash_value) {
        num_hash++;
    }

    // Used for second level match
    //   so only need location and hash value
    oridata[i].sig4.reserve(sig3.size());
    for (int j = 0; j < sig3.size(); j++) {
        const Signature_inter& s = sig3[j];
        oridata[i].sig4.push_back({s.location, s.hash_value});
    }

    sig3.resize(num_hash+1);
    oridata[i].num_sampled_hash = num_hash+1;
    return sig3;
}

int get_top_k_vec(const std::vector<int>& hash_bucket, double topk = TOPK) {
    std::vector<int>  num_element_per_hash;
    num_element_per_hash.reserve(hash_bucket.size()-1);
    for (int i = 1; i < hash_bucket.size(); i++) {
        num_element_per_hash.push_back(hash_bucket[i] - hash_bucket[i-1]);
    }

    sort(num_element_per_hash.begin(), num_element_per_hash.end());
    int index = (int) ((1-topk) * num_element_per_hash.size());
    return num_element_per_hash[index];
}
