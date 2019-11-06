#pragma once

#include "util.h"
#include <utility>
#include "build_index.h"

std::vector<Signature_inter> generate_signature_vec() {
    // Controls the generation of smoothq grams
    srand(RANDOM_SEED);
    // Controls the smq hash scheme
    int seed = HASH_SEED;

    // Build randomness
    std::vector<int> hash_lsh = generate_random_sampling_positions(SIZE_CGK, SIZE_SMQ);
    int **p = generate_random_binary(SIZE_CGK);

    double current_cpu = get_cpu_time();
    double current_time = time(NULL);
    std::vector<Signature_inter> ht;

    #pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS)
    for (int i = 0; i < oridata.size(); i++) {
        std::vector<Signature_inter> sig3 = generate_ith_signature(i, p, hash_lsh, seed);

        // Only allow single thread to operate at one time
        #pragma omp critical
        ht.insert(ht.end(), sig3.begin(), sig3.end());

        sig3.clear();
    }
    fprintf(stderr, "Finished generating signatures in %fs cpu time \n", get_cpu_time() - current_cpu);
    fprintf(stderr, "                               in %fs real time \n", time(NULL) - current_time);
    return ht;
}

std::vector<Signature_inter> generate_table_vec(std::vector<Signature_inter>& ht) {
    // Sort initial ht vector
    double current_cpu = get_cpu_time();
    std::sort(ht.begin(), ht.end(), sortSiginter);
    unsigned long ht_size = ht.size();

    std::vector<int> hash_bucket;
    hash_bucket.push_back(0);
    HASH_TYPE current_hash = ht[0].hash_value;
    for (int i = 0; i < ht_size; i++) {
        if (ht[i].hash_value != current_hash) {
            current_hash = ht[i].hash_value;
            hash_bucket.push_back(i);
        }
    }
    hash_bucket.push_back(ht_size);
    fprintf(stderr, "Initial vec hash table contains %lu hash values actually \n", hash_bucket.size()-1);
    fprintf(stderr, "                       contains %lu elements actually \n", ht_size);

    // filter initial hash table into htable
    std::vector<Signature_inter> htable;
    htable.reserve(ht.size());
    FREQUENT = get_top_k_vec(hash_bucket);
    fprintf(stderr, "Frequency threshold is set to be %d with top %f \n",
            FREQUENT, TOPK);

    unsigned long num_singular_hash = 0;
    unsigned long num_singular_elements = 0;
    unsigned long num_frequent_hash = 0;
    unsigned long num_frequent_elements = 0;
    for (int i = 0; i < hash_bucket.size()-1; i++) {
        int num_signature_same_hash = hash_bucket[i+1] - hash_bucket[i];
        // filtering step
        if (num_signature_same_hash <= SINGULAR) {
            num_singular_hash++;
            num_singular_elements += num_signature_same_hash;
        }
        else if (num_signature_same_hash > FREQUENT) {
            num_frequent_hash++;
            num_frequent_elements += num_signature_same_hash;
        }
        else htable.insert(htable.end(), ht.begin()+hash_bucket[i], ht.begin()+hash_bucket[i+1]);
    }

    /*
    fprintf(stderr, "Filtered vec hash table contains %lu hash values, left with %f%% of initial\n",
            hash_bucket.size()-1-num_singular_hash-num_frequent_hash,
            100*(1-(double) (num_singular_hash+num_frequent_hash) / (hash_bucket.size()-1)));
    fprintf(stderr, "  Filtered out %lu singular hash values, %f%% of initial\n",
            num_singular_hash, 100*(double)num_singular_hash/(hash_bucket.size()-1));
    fprintf(stderr, "  Filtered out %lu frequent hash values, %f%% of initial\n",
            num_frequent_hash, 100*(double)num_frequent_hash/(hash_bucket.size()-1));
    */
    fprintf(stderr, "Filtered vec hash table contains %lu elements, left with %f%% of elements \n",
            htable.size(), 100* (double) htable.size() / ht_size);
    fprintf(stderr, "  Filtered out %lu singular elements, %f%% of initial\n",
            num_singular_elements, 100*(double)num_singular_elements/ht_size);
    fprintf(stderr, "  Filtered out %lu frequent elements, %f%% of initial\n",
            num_frequent_elements, 100*(double)num_frequent_elements/ht_size);
    fprintf(stderr, "Finished building final hash table in %f s cpu \n", get_cpu_time()-current_cpu);
    ht.clear();

    return htable;
}
