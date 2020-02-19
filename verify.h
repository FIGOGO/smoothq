#pragma once

#include "util.h"
#include "density.h"

#include <iterator>

// estimate overlap length by the median x-y pair
int estimate_overlap_length(const std::vector<Match>& valid_matches, const int i, const int j, const int offset_median) {
    int x = 0, y = 0;
    int x_end = oridata[i].len, y_end = oridata[j].len;

    for (int k = 0; k < valid_matches.size(); k++) {
        int offset = valid_matches[k].loc1 - valid_matches[k].loc2;
        if (offset == offset_median) {
            x = valid_matches[k].loc1;
            y = valid_matches[k].loc2;
            break;
        }
    }
    int ovl_before = std::min(x, y);
    int ovl_after = std::min(x_end-x, y_end-y);
    int estimated_ovl = ovl_before + ovl_after;
    if (DEBUG) fprintf(stderr, "Estimated overlap length is %d, shift may be %d for reads with length %d and %d, "
                               "median offset is %d\n",
                       estimated_ovl, (int) (estimated_ovl*DATA_ERROR_RATE), x_end, y_end, offset_median);
    return estimated_ovl;
}

bool first_phase_check(int i, int j, std::vector<Match> &matches, int *offset_median, int shift_size = SHIFT,
                       int window_size = WSIZE) {
    double current_cpu_time = omp_get_wtime();
    std::vector<Match> validmatches;
    validmatches.reserve(matches.size());

    if (DEBUG) {
        fprintf(stderr, "\nWorking on data %d and %d\n", oridata[i].index, oridata[j].index);
        fprintf(stderr, "find offset 1st\n");
    }
    *offset_median = find_offset(matches, validmatches, shift_size);
    // If the matches after offset are sufficient
    if (validmatches.size() > 20*THRESHOLD)
        return true;
    // If the matches after offset are not enough
    if (validmatches.size() < THRESHOLD)
        return false;

    if (DEBUG) fprintf(stderr, "find overlap 1st\n");
    int validmatches_range = find_overlap_area(validmatches, window_size);
    if (EFFICIENCY_CHECK) {
        double now_cpu_time = omp_get_wtime();
#pragma omp critical
        TIME_1st_CHECK += (now_cpu_time-current_cpu_time);
    }
    return (validmatches_range >= MIN_SPAN);
//    return find_match_in_wsize(validmatches, window_size);
}

std::vector<Match> calmatches(int i, int j, int offset_median, int estimated_overlap)
{
    const std::vector<Signature_intra>& x = oridata[i].sig4;
    const std::vector<Signature_intra>& y = oridata[j].sig4;
    if (oridata[i].hash_bucket.empty()) {
        oridata[i].hash_bucket = get_hash_bucket(x);
    }
    if (oridata[j].hash_bucket.empty()) {
        #pragma omp critical
        oridata[j].hash_bucket = get_hash_bucket(y);
    }
    const std::vector<HashBucket>& bi = oridata[i].hash_bucket;
    const std::vector<HashBucket>& bj = oridata[j].hash_bucket;
    std::vector<Match> output;
    output.reserve(x.size());

    // calculate offset range allowance
    int min_offset = offset_median - std::max((int) (estimated_overlap*DATA_ERROR_RATE), SHIFT_final);
    int max_offset = offset_median + std::max((int) (estimated_overlap*DATA_ERROR_RATE), SHIFT_final);

    // init counters
    int i1 = 0;
    int i2 = 0;

    // perform merge operation to get the shift and the kmer count
    while (i1 < bi.size() && i2 < bj.size()) {

        // get the values in the array
        HASH_TYPE hash1 = bi[i1].hash_value;
        HASH_TYPE hash2 = bj[i2].hash_value;

        if (hash1 < hash2)
            i1++;
        else if (hash2 < hash1)
            i2++;
        else {
            // Try to capture tandem repeats
            if (bi[i1].end-bi[i1].begin > MAX_NUM_PER_SIG && bj[i2].end-bj[i2].begin > MAX_NUM_PER_SIG) {
                i1++;
                i2++;
                continue;
            }
            for (int k = bi[i1].begin; k < bi[i1].end; k++) {
                for (int l = bj[i2].begin; l < bj[i2].end; l++) {
                    int offset = x[k].location - y[l].location;

                    if (offset < max_offset && offset > min_offset) {
                        output.push_back({ x[k].location, y[l].location });
                    }
                }
            }
            i1++;
            i2++;
        }
    }
    if(DEBUG) {
//        print_match_info(output);
        fprintf(stderr, "Total number of matched pairs: %d\n", output.size());
    }
    return output;
}

// Add flexible shift_final based on estimated overlap length
void second_phase_match(int i, int j, int offset_median, int estimated_overlap) {
    double current_cpu_time = omp_get_wtime();
    double cpu_time = omp_get_wtime();

    std::vector<Match> newmatches, newvalidmatches;
    newmatches = calmatches(i, j, offset_median, estimated_overlap);
    if (newmatches.size() <= THRESHOLD) return;

    if (EFFICIENCY_CHECK) {
        double now_cpu_time = omp_get_wtime();
        #pragma omp critical
        TIME_CAL_MATCH += (now_cpu_time-cpu_time);
    }

    if(DEBUG) fprintf(stderr, "find offset 2nd\n");

    find_offset(newmatches, newvalidmatches, std::max((int) (estimated_overlap*DATA_ERROR_RATE_r2), SHIFT_final));

    if (NUMBER_CHECK) {
        #pragma omp critical
        NUM_INTO_EXTEND++;
    }

//    if (newvalidmatches.size() > (int) (std::min(oridata[i].len, oridata[j].len)*0.5)) {
//        find_offset(newvalidmatches, newvalidmatches2, std::min(400, (int) (0.8 * std::max((int) (estimated_overlap*DATA_ERROR_RATE), SHIFT_final))));
//    }

    if (EFFICIENCY_CHECK) {
        double now_cpu_time = omp_get_wtime();
        #pragma omp critical
        TIME_2nd_CHECK += (now_cpu_time-current_cpu_time);
    }

    current_cpu_time = omp_get_wtime();
    extendmatch(i, j, newvalidmatches);

    if (EFFICIENCY_CHECK) {
        double now_cpu_time = omp_get_wtime();
        #pragma omp critical
        TIME_EXTEND += (now_cpu_time-current_cpu_time);
    }
}


int verify(int i, int j, std::vector<Match> &matches)
{
//    fprintf(stderr, "First phase check for %d and %d\n", oridata[i].index, oridata[j].index);
    int offset_median = 0;
    // If passed first phase check, enter second phase match
    if (first_phase_check(i, j, matches, &offset_median)) {
        if (NUMBER_CHECK) {
            #pragma omp critical
            NUM_INTO_2nd_CHECK++;
        }
        int estimated_overlap = estimate_overlap_length(matches, i, j, offset_median);
//        fprintf(stderr, "Second phase check for %d and %d\n", oridata[i].index, oridata[j].index);
        second_phase_match(i, j, offset_median, estimated_overlap);
    }

    return 0;
}

void join(std::vector<Signature_inter>& htable) {
    std::unordered_map<HASH_TYPE, std::pair<int, int>> hash_index = generate_hash_index(htable);

    int num_verify_threads = NUM_THREADS;
    if (DEBUG) num_verify_threads = 1;

    #pragma omp parallel for schedule(dynamic) num_threads(num_verify_threads)
    for (int i = 0; i < NUM_READS_DATA1; i++) {
        if (oridata[i].fwrv == 0) {
            double current_cpu_time = omp_get_wtime();
            std::vector<RawMatch> raw_matches;
            for (int j = 0; j < oridata[i].num_sampled_hash; j++) {
                const Signature_intra& s1 = oridata[i].sig4[j];
                if (hash_index.find(s1.hash_value) == hash_index.end()) {
                    continue;
                }
                const std::pair<int, int>& index_pair = hash_index[s1.hash_value];
                Signature_inter base_sig = {i+2, s1.location, s1.hash_value};
                std::vector<Signature_inter>::iterator begin_iterator = std::lower_bound(
                        htable.begin() + index_pair.first, htable.begin() + index_pair.second, base_sig, sortSiginter);
                int begin = begin_iterator - htable.begin();
                for (int k = begin; k < index_pair.second; k++) {
                    const Signature_inter& s2 = htable[k];
                    if (NUM_READS_DATA2 != 0 && s2.id < NUM_READS_DATA1) continue;
                    if (edit_early_quit(oridata[i].cstr+s1.location, oridata[s2.id].cstr+s2.location, SIZE_Q, ED_THRESHOLD) != -1) {
                        raw_matches.push_back({s2.id, s1.location, s2.location});
                    }
                }
            }
            sort(raw_matches.begin(), raw_matches.end(), sortRawMatch);
            std::vector<MatchBucket> match_bucket = get_match_index(raw_matches);

            if (NUMBER_CHECK) {
                #pragma omp critical
                NUM_INTO_VERIFY += match_bucket.size();
            }

            if (EFFICIENCY_CHECK) {
                double now_cpu_time = omp_get_wtime();
                #pragma omp critical
                TIME_BEFORE_VERIFY += (now_cpu_time-current_cpu_time);
            }

            current_cpu_time = omp_get_wtime();
            for (auto bucket : match_bucket) {
                std::vector<Match> matches;
                matches.reserve(bucket.end - bucket.begin);
                for (int k = bucket.begin; k < bucket.end; k++) {
                    matches.push_back({raw_matches[k].loc1, raw_matches[k].loc2});
                }
                verify(i, bucket.id, matches);
                if (NUMBER_CHECK) {
                    #pragma omp critical
                    NUM_INTO_VERIFY++;
                }
            }

            if (EFFICIENCY_CHECK) {
                double now_cpu_time = omp_get_wtime();
                #pragma omp critical
                TIME_IN_VERIFY += (now_cpu_time-current_cpu_time);
            }
        }
    }
}
