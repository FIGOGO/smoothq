#pragma once

bool step_qualification(const Match& x, const Match& x_neighbor,
                        const double thr1 = EXTEND_THRESHOLD) {
    int shift_loc1 = x.loc1 - x_neighbor.loc1;
    int shift_loc2 = x.loc2 - x_neighbor.loc2;
    // Check crossing before extending
    if (shift_loc1 * shift_loc2 <= 0) return false;

    int step = (abs(shift_loc1) + abs(shift_loc2)) / 2;
    if (step >= EXTENDSTEP) return false;

    int shift_diff = abs(shift_loc1 - shift_loc2);
    int shift_max = std::max(abs(shift_loc1), abs(shift_loc2));
    double shift_percent = (double) shift_diff / shift_max;

    // Qualification to connect
    if (step < EXTENDSTEP && shift_diff < 20) return true;
    // step in range [0, extendstep)
    if (step < EXTENDSTEP && shift_percent < thr1) return true;
    return false;
}

void print_extend_info(const Match& x, const Match& x_neighbor) {
    int shift_loc1 = x.loc1 - x_neighbor.loc1;
    int shift_loc2 = x.loc2 - x_neighbor.loc2;
    int step = (abs(shift_loc1) + abs(shift_loc2)) / 2;
    int shift_diff = abs(shift_loc1 - shift_loc2);
    int shift_max = std::max(abs(shift_loc1), abs(shift_loc2));
    double shift_percent = (double) shift_diff / shift_max;

    fprintf(stderr, "  x: [%d, %d], x_neighbor: [%d, %d], "
                    "shift_loc1 %d, shift_loc2 %d, shift_diff %d, "
                    "step %d, shift_percent %.3f \n",
            x.loc1, x.loc2, x_neighbor.loc1, x_neighbor.loc2,
            shift_loc1, shift_loc2, shift_diff,
            step, shift_percent);
}

bool extend_one_step(Match& x, const Match& x_neighbor) {
    if (step_qualification(x, x_neighbor)) {
        x = x_neighbor;
        return true;
    } else {
        if (DEBUG) print_extend_info(x, x_neighbor);
        return false;
    }
}
