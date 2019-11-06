#pragma once

#include "util.h"
#include "extend_match.h"

void find_dense_window(const std::vector<int>& data, const int wsize, int& begin, int& end) {
    const int n = data.size();

    std::vector<int> position(n, 0);

    begin = 0;
    end = 0;
    int maximum_density = 1;
    for (int i = 1; i < n; i++) {
        // Fill in position[i]
        int j = 0;
        for (j = position[i-1]; j < i; j++) {
            if (data[i] - data[j] <= wsize) {
                position[i] = j;
                break;
            }
        }
        if (i == j) position[i] = i;

        // Find max density
        int density = i - position[i] + 1;
        if (density > maximum_density) {
            begin = position[i];
            end = i;
            maximum_density = density;
        }
    }

    if (DEBUG) {
        fprintf(stderr, "Out of %d chosen %d data points in index [%d, %d]."
                        " Begins at %d, ends at %d, actural range %d \n",
                        data.size(), end-begin+1, begin, end,
                        data[begin], data[end], data[end] - data[begin]);

    }
}

int find_offset(std::vector<Match>& matches, std::vector<Match>& valid_matches, int wsize) {
    std::vector<int> offsets(matches.size(), 0);
    for (int i = 0; i < matches.size(); i++) {
        offsets[i] = matches[i].loc1 - matches[i].loc2;
    }
    sort(offsets.begin(), offsets.end());

    int begin, end;
    find_dense_window(offsets, wsize, begin, end);

    int low_offset = offsets[begin];
    int high_offset = offsets[end];
    for (int i = 0; i < matches.size(); i++) {
        int cur_offset = matches[i].loc1 - matches[i].loc2;
        if (cur_offset >= low_offset && cur_offset <= high_offset) {
            valid_matches.push_back(matches[i]);
        }
    }

    // return median inside the window size
    return offsets[(begin+end) / 2];
}

int find_overlap_area(const std::vector<Match>& valid_matches, int wsize) {
    std::vector<int> locations(valid_matches.size(), 0);
    for (int i = 0; i < valid_matches.size(); i++) {
        locations[i] = valid_matches[i].loc1;
    }
    sort(locations.begin(), locations.end());

    int begin, end;
    find_dense_window(locations, wsize, begin, end);

    int low_location = locations[begin];
    int high_location = locations[end];

    return high_location - low_location + 1;
}

void print_match_info(std::vector<Match>& matches) {
    sort(matches.begin(), matches.end(), sortMatch);
    for (int i = 0; i < matches.size(); i++) {
        fprintf(stderr, "[%d %d]_%d, ",
                matches[i].loc1, matches[i].loc2,
                matches[i].loc1 - matches[i].loc2);
    }
    fprintf(stderr, "\n");
}

std::vector<HashBucket> get_hash_bucket(const std::vector<Signature_intra>& sig) {
    std::vector<HashBucket> bucket;
    bucket.reserve(sig.size());
    HASH_TYPE current_hash = sig[0].hash_value;
    int index = 0;
    for (int i = 0; i < sig.size(); i++) {
        if (sig[i].hash_value != current_hash) {
            bucket.push_back({current_hash, index, i});
            index = i;
            current_hash = sig[i].hash_value;
        }
    }

    // Add the last hash value
    bucket.push_back({current_hash, index, (int) sig.size()});
    return bucket;
}


std::unordered_map<HASH_TYPE, std::pair<int, int>> generate_hash_index(std::vector<Signature_inter>& htable) {
    if (htable.empty()) return {};

    // Build hash_value to index pair table
    std::unordered_map<HASH_TYPE, std::pair<int, int>> hash_index;
    int index_begin = 0;
    HASH_TYPE current_hash = htable[0].hash_value;
    for (int i = 0; i < htable.size(); i++) {
        if (current_hash != htable[i].hash_value) {
            hash_index.insert({current_hash, std::make_pair(index_begin, i)});
            index_begin = i;
            current_hash = htable[i].hash_value;
        }
    }
    hash_index.insert({current_hash, std::make_pair(index_begin, htable.size())});
    return hash_index;
}

std::vector<MatchBucket> get_match_index(std::vector<RawMatch>& raw_matches) {
    if (raw_matches.size() == 0) return {};

    std::vector<MatchBucket> match_bucket;
    match_bucket.reserve(raw_matches.size());

    int current_id = raw_matches[0].id;
    int index_begin = 0;
    for (int i = 1; i < raw_matches.size(); i++) {
        if (current_id != raw_matches[i].id) {
            if (i - index_begin >= THRESHOLD) {
                match_bucket.push_back({current_id, index_begin, i});
            }
            index_begin = i;
            current_id = raw_matches[i].id;
        }
    }

    match_bucket.push_back({current_id, index_begin, (int) raw_matches.size()});
    return match_bucket;
}