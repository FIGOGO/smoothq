#pragma once

#include "util.h"
#include "step_extension.h"

void print_window_merge_info(const Window& w, const Window& x) {
    int shift_loc1 = x.r1_loc1 - w.r1_loc2;
    int shift_loc2 = x.r2_loc1 - w.r2_loc2;
    int step = (abs(shift_loc1) + abs(shift_loc2)) / 2;
    int shift_diff = abs(shift_loc1 - shift_loc2);
    int shift_max = std::max(abs(shift_loc1), abs(shift_loc2));
    double shift_percent = (double) shift_diff / shift_max;

    fprintf(stderr, "  w1: [%d, %d] -> [%d, %d], w2: [%d, %d] -> [%d, %d], "
                    "shift_loc1 %d, shift_loc2 %d, shift_diff %d, "
                    "step %d, shift_percent %.3f \n",
            w.r1_loc1, w.r2_loc1, w.r1_loc2, w.r2_loc2,
            x.r1_loc1, x.r2_loc1, x.r1_loc2, x.r2_loc2,
            shift_loc1, shift_loc2, shift_diff,
            step, shift_percent);
}

// Only care the step size; shift_percent now considered
bool windows_merge_qualification(const Window& w, const Window& x){
    int max_window_gap = MAX_WINDOW_GAP;

    int shift_loc1 = x.r1_loc1 - w.r1_loc2;
    int shift_loc2 = x.r2_loc1 - w.r2_loc2;
//    if (shift_loc1 * shift_loc2 <= 0) return false;
    if (x.r1_loc2 < w.r1_loc1 || x.r2_loc2 < w.r2_loc1) return false;

    int step = (abs(shift_loc1) + abs(shift_loc2)) / 2;
    if (step > max_window_gap) return false;

    int large_window_size = std::max(w.size(), x.size());
    if (step < large_window_size) return true;

    int shift_diff = abs(shift_loc1 - shift_loc2);
    int shift_max = std::max(abs(shift_loc1), abs(shift_loc2));
    double shift_percent = (double) shift_diff / shift_max;

    return (shift_percent < MAX_WINDOW_PERCENT && step < 1.5*large_window_size);
}

// Only try to merge adjacent windows
std::vector <Window> merge_adjacent_windows(std::vector<Window>& windows) {
    if (windows.size() <= 1) return windows;

    sort(windows.begin(), windows.end(), sortWindow);
    std::vector<Window> new_windows;
    Window begin_window = windows[0];

    for (int i = 1; i < windows.size(); i++) {
        if (DEBUG) print_window_merge_info(begin_window, windows[i]);

        if (windows_merge_qualification(begin_window, windows[i])){
            begin_window.r1_loc2 = windows[i].r1_loc2;
            begin_window.r2_loc2 = windows[i].r2_loc2;
        }
        else {
            new_windows.push_back(begin_window);
            begin_window = windows[i];
        }
    }
    new_windows.push_back(begin_window);

    if (DEBUG) {
        fprintf(stderr, "Windows after merge: \n");
        for (int i = 0; i < new_windows.size(); i++) {
            const Window& w = new_windows[i];
            fprintf(stderr, "  [%d %d] -> [%d %d]. Offset: %d, num sig: %d, size: %d\n",
                    w.r1_loc1, w.r2_loc1, w.r1_loc2, w.r2_loc2,
                    w.r1_loc1 - w.r2_loc1, w.num_sig, w.size());
        }
    }
    return new_windows;
}

// Return after two rounds of merge_adjacent_windows
Window* merge_overlap_windows(std::vector<Window>& windows) {
    if (windows.size() == 0) return nullptr;
    if (windows.size() == 1) return &(windows[0]);

    int window_size = windows.size();
    if (DEBUG) fprintf(stderr, "To merge adjacent windows 1st round:\n");
    windows = merge_adjacent_windows(windows);

    if (windows.size() > 1 and windows.size() < window_size) {
        if (DEBUG) fprintf(stderr, "To merge adjacent windows 2nd round:\n");
        windows = merge_adjacent_windows(windows);
    }

    // Find the window with largest overlap area
    int max_ovl_length = 0;
    int max_window_index = -1;
    Window* res = new Window();
    for (int i = 0; i < windows.size(); i++) {
        Window& w = windows[i];
        int ovl_length = (w.r1_loc2-w.r1_loc1 + w.r2_loc2-w.r2_loc1) / 2;
        if (ovl_length > max_ovl_length) {
            max_ovl_length = ovl_length;
            max_window_index = i;
            *res = w;
        }
    }

    // Possible extend to left
    for (int i = max_window_index-1; i >= 0; i--) {
        if (windows_merge_qualification(windows[i], *res)){
            res->r1_loc1 = windows[i].r1_loc1;
            res->r2_loc1 = windows[i].r2_loc1;
        }
    }

    // Possible extend to left
    for (int i = max_window_index+1; i < windows.size(); i++) {
        if (windows_merge_qualification(*res, windows[i])){
            res->r1_loc2 = windows[i].r1_loc2;
            res->r2_loc2 = windows[i].r2_loc2;
        }
    }

    if (DEBUG){
        fprintf(stderr, "Merged window should be");
        fprintf(stderr, " [%d %d] -> [%d %d]\n", res->r1_loc1, res->r2_loc1, res->r1_loc2, res->r2_loc2);
    }
    return res;
}

void merge_matches_to_windows (std::vector<Match>& matches, std::vector<Window>& windows, int max_num_failure = 5) {
    Match xbegin = matches[0];

    int num_failure = 0;
    int x_begin_loc1 = xbegin.loc1;
    int x_begin_loc2 = xbegin.loc2;
    int i = 1;
    int num_sig = 0;
    if (DEBUG) {
        fprintf(stderr, "Begin with match [%d, %d]\n", xbegin.loc1, xbegin.loc2);
    }
    std::vector<Match> window_path;
    window_path.push_back(xbegin);

    Window w;
    w.add_offset(xbegin);
    while (i < matches.size())
    {
        if (extend_one_step(xbegin, matches[i])) {
            num_failure = 0;
            num_sig++;
            w.add_offset(matches[i]);
            if (DEBUG) {
                window_path.push_back(matches[i]);
            }
        }
        else if (num_failure < max_num_failure-1) {
            num_failure++;
        }
        else {
            num_failure = 0;

            w.r1_loc1 = x_begin_loc1;
            w.r2_loc1 = x_begin_loc2;
            w.r1_loc2 = xbegin.loc1;
            w.r2_loc2 = xbegin.loc2;
            w.num_sig = num_sig;
            if (abs(w.r1_loc1 - w.r1_loc2) > MIN_WINDOW_SIZE-SIZE_Q && w.sig_density() > MIN_WINDOW_DENSITY && w.num_sig > 2)
                windows.push_back(w);

            num_sig = 0;
            i = i-max_num_failure+1;
            xbegin = matches[i];

            if (DEBUG) {
                fprintf(stderr, "  Window path: ");
                for (int j = 0; j < window_path.size(); j++) {
                    const Match& m = window_path[j];
                    fprintf(stderr,"[%d %d], ", m.loc1, m.loc2);
                }
                fprintf(stderr, "\n");
                fprintf(stderr, "  Num_sig: %d, num offsets: %d, size: %d, density %f \n", w.num_sig, w.offsets.size(), w.size(), w.sig_density());
                window_path.clear();
                window_path.push_back(xbegin);

                fprintf(stderr, "Begin with match [%d, %d]\n", xbegin.loc1, xbegin.loc2);
            }

            w.clear();
            w.add_offset(xbegin);

            x_begin_loc1 = xbegin.loc1;
            x_begin_loc2 = xbegin.loc2;
        }
        i++;
    }

    // Insert the last window
    w.r1_loc1 = x_begin_loc1;
    w.r2_loc1 = x_begin_loc2;
    w.r1_loc2 = xbegin.loc1;
    w.r2_loc2 = xbegin.loc2;
    w.num_sig = num_sig;
    if (abs(w.r1_loc1 - w.r1_loc2) > MIN_WINDOW_SIZE-SIZE_Q && w.sig_density() > MIN_WINDOW_DENSITY && w.num_sig > 2)
        windows.push_back(w);
    if (DEBUG) {
        fprintf(stderr, "  Window path: ");
        for (int j = 0; j < window_path.size(); j++) {
            const Match& m = window_path[j];
            fprintf(stderr,"[%d %d], ", m.loc1, m.loc2);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "  Num_sig: %d, size: %d, density %f \n", w.num_sig, w.size(), w.sig_density());
        window_path.clear();
    }
}
