#pragma once

#include "util.h"
#include "window_extention.h"

void output_overlap(int i, int j, Match& xbegin, Match& xend, double density) {
    int a1, a2, b1, b2;
    a1 = xbegin.loc1;
    a2 = xend.loc1;
    b1 = xbegin.loc2;
    b2 = xend.loc2;

    a2 = std::min((int)oridata[i].len-1, a2 + SIZE_Q);
    b2 = std::min((int)oridata[j].len-1, b2 + SIZE_Q);

    if (oridata[j].fwrv)
    {
        int tmp = b2;
        b2 = oridata[j].len - b1;
        b1 = oridata[j].len - tmp;
    }


    if (a2 - a1 > 200 || b2 - b1 > 200) {
        if (DEBUG) fprintf(stderr, "output: \n");
        if (EFFICIENCY_CHECK) {
            #pragma omp critical
            NUM_INTO_OUTPUT++;
        }
        #pragma omp critical
        std::cout << oridata[i].index
                  << " " << oridata[j].index
                  << " " << oridata[i].tag << "_" << oridata[j].tag
                  << " " << density
                  << " " << oridata[i].fwrv
                  << " " << a1 << " " << a2 << " " << oridata[i].len
                  << " " << oridata[j].fwrv
                  << " " << b1 << " " << b2 << " " << oridata[j].len
                  << std::endl;
    }
}

// Extending through match -> window -> segment
void extendmatch(int index1, int index2, std::vector<Match>& matches) {
    sort(matches.begin(), matches.end(), sortMatch);
    std::vector<Window> windows;
    merge_matches_to_windows(matches, windows);

    if(DEBUG) {
        fprintf(stderr, "merge matches to windows for read %d and %d\n", oridata[index1].index, oridata[index2].index);
        fprintf(stderr, "There are %d windows created\n",windows.size());
        for (int i = 0; i < windows.size(); i++) {
            Window w = windows[i];
            fprintf(stderr, "  [%d %d] -> [%d %d]", w.r1_loc1, w.r2_loc1, w.r1_loc2, w.r2_loc2);
            fprintf(stderr, "    Num_sig: %d, size: %d, density %f \n", w.num_sig, w.size(), w.sig_density());
        }
    }

    if (windows.size() == 0) return;

    Window* w;
    if (windows.size() == 1) w = &windows[0];
    else w = merge_overlap_windows(windows);

    if (w == nullptr) return;

    if (DEBUG) {
        fprintf(stderr, "merge overlap windows for read %d and %d\n", oridata[index1].index, oridata[index2].index);
        fprintf(stderr, "Merged window [%d %d] -> [%d %d]", w->r1_loc1, w->r2_loc1, w->r1_loc2, w->r2_loc2);
        fprintf(stderr, "    Num_sig: %d, size: %d, density %f \n", w->num_sig, w->size(), w->sig_density());
    }


    // Current xbegin and xend
    Match xbegin;
    xbegin.loc1 = w->r1_loc1;
    xbegin.loc2 = w->r2_loc1;
    Match xend;
    xend.loc1 = w->r1_loc2;
    xend.loc2 = w->r2_loc2;

    output_overlap(index1, index2, xbegin, xend, w->sig_density());
}
