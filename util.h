#pragma once

#include <unordered_map>
#include <tuple>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

typedef unsigned long long HASH_TYPE;

// Useful data structure
struct Signature_inter
{
    int id;
    int location;
    HASH_TYPE hash_value;
};

struct Signature_intra
{
    int location;
    HASH_TYPE hash_value;
};

struct Match
{
    int loc1;
    int loc2;
};

struct HashBucket {  // used in recalculating matches
    HASH_TYPE hash_value;
    int begin;
    int end;
};

struct RawMatch {
    int id;
    int loc1;
    int loc2;
};

struct MatchBucket {
    int id;
    int begin;
    int end;
};

struct Data {
	int index;
    std::string str; //initialize in readdata()
    const char * cstr;
    unsigned long len;
    std::string tag;
    int num_sampled_hash;

    int fwrv; // fwd=0 rev=1;
    std::vector<Signature_intra> sig4;
    std::vector<HashBucket> hash_bucket;
};

class Window {
public:
    int r1_loc1;
    int r1_loc2;
    int r2_loc1;
    int r2_loc2;
    int num_sig;
    std::vector<int> offsets;

    Window(){
        this->num_sig = 0;
    }

    // Copy constructor
    Window(const Window& w) {
        this->r1_loc1 = w.r1_loc1;
        this->r2_loc1 = w.r2_loc1;
        this->r1_loc2 = w.r1_loc2;
        this->r2_loc2 = w.r2_loc2;
        this->num_sig = w.num_sig;
        this->offsets = w.offsets;
    }

    int size() const { return (this->r1_loc2-this->r1_loc1 + this->r2_loc2-this->r2_loc1) / 2; }

    double sig_density() { return ((double) this->num_sig) / this->size(); }

    // Extend the window
    void extend(const Window& w) {
        this->r1_loc2 = w.r1_loc2;
        this->r2_loc2 = w.r2_loc2;
        this->num_sig += w.num_sig;
    }

    void extend_match(const Match& m) {
        if (m.loc1 < this->r1_loc1) extend_match_ahead(m);
        if (m.loc1 > this->r1_loc2) extend_match_behind(m);
    }

    void extend_match_ahead(const Match& m) {
        this->r1_loc1 = m.loc1;
        this->r2_loc1 = m.loc2;
        this->num_sig++;
    }

    void extend_match_behind(const Match& m) {
        this->r1_loc2 = m.loc1;
        this->r2_loc2 = m.loc2;
        this->num_sig++;
    }

    int add_offset(const Match& m) {
        int offset = m.loc1 - m.loc2;
        this-> offsets.push_back(offset);
        return offset;
    }

    void clear() {
        this->r1_loc1 = 0;
        this->r2_loc1 = 0;
        this->r1_loc2 = 0;
        this->r2_loc2 = 0;
        this->num_sig = 0;
        this->offsets.clear();
    }
};

// Helper functions
double get_cpu_time(){
	return (double)clock() / CLOCKS_PER_SEC;
}

inline bool sortSiginter(const Signature_inter &a, const Signature_inter &b)
{
	if (a.hash_value == b.hash_value) return a.id < b.id;
	else return a.hash_value < b.hash_value;
}

inline bool sortMatch(const Match &a, const Match &b)
{
	// smallest comes first
	return a.loc1 < b.loc1;
}

inline bool sortWindow(const Window &x, const Window &y) {
    return x.r1_loc1 < y.r1_loc1;
}


inline bool sortRawMatch (const RawMatch& a, const RawMatch& b)
{ return a.id < b.id; }



// Input data
std::vector<Data> oridata; //input data

std::unordered_map<char, char> trans({ { 'A', 'T' }, { 'T', 'A' }, { 'G', 'C' }, { 'C', 'G' } });
std::unordered_map<char, int> dict({ { 'A', 0 }, { 'C', 1 }, { 'G', 2 }, { 'T', 3 }, { 'N', 4 } });

//Parameters to be tuned
int RANDOM_SEED = 100;
int HASH_SEED = 100;

// Parameters for input file
const int MIN_LEN = 500; // minimum length of input to consider

// Parameters for CGK embedding
int SIZE_Q = 14; //len of qgram
int SIZE_CGK = 35;
int SIZE_SMQ = 16;

double SAMPLING_RATE = 0.2;// signature selection rate
double TOPK = 0.005; //frequency filtering threshold
int FREQUENT = 5e3;
int SINGULAR = 1;
int NUM_THREADS = 16;
int DEBUG = 0;

//parameters for verification
// First phase checking
int THRESHOLD = 5;
int ED_THRESHOLD = 2;
const int SHIFT = 200;
const int WSIZE = 500;
const int MIN_SPAN = 3;

// Second phase checking
const int SHIFT_final = 500;
const double DATA_ERROR_RATE = 0.15;
const double DATA_ERROR_RATE_r2 = 0.1;
const int MAX_NUM_PER_SIG = 5;  // if one signature appears at too many locations, maybe tandem repeats

// Step extension
const int    EXTENDSTEP       = 200; // maximum step to extend in verification
const double EXTEND_THRESHOLD = 0.15;

// Window extending
const int MIN_WINDOW_SIZE = 20;
const int MAX_WINDOW_GAP = 2000;
const double MIN_WINDOW_DENSITY = 0.00;
const double MAX_WINDOW_PERCENT = 0.4;

// Efficiency check
const int EFFICIENCY_CHECK = 0;
const int NUMBER_CHECK = 0;
unsigned long long NUM_INTO_VERIFY = 0;
unsigned long long NUM_INTO_2nd_CHECK = 0;
unsigned long long NUM_INTO_EXTEND = 0;
unsigned long long NUM_INTO_OUTPUT = 0;
double TIME_BEFORE_VERIFY = 0;
double TIME_IN_VERIFY = 0;
double TIME_1st_CHECK = 0;
double TIME_CAL_MATCH = 0;
double TIME_2nd_CHECK = 0;
double TIME_EXTEND = 0;

// Output format
std::string OUT_FORMAT = "m4";

void print_working_env() {
    fprintf(stderr, "Parameter list: \n");
    fprintf(stderr, "  Q-gram size:         %d \n", SIZE_Q);
    fprintf(stderr, "  CGK-embed size:      %d \n", SIZE_CGK);
    fprintf(stderr, "  Smooth-q-gram size:  %d \n", SIZE_SMQ);
    fprintf(stderr, "  Sampling rate:       %f \n", SAMPLING_RATE);
    fprintf(stderr, "  Filter:              %f \n", TOPK);
    fprintf(stderr, "  Min matches size:    %d \n", THRESHOLD);
    fprintf(stderr, "  Max edit distance:   %d \n", ED_THRESHOLD);
    fprintf(stderr, "  Num threads:         %d \n", NUM_THREADS);
}


int number_check() {
    if(NUMBER_CHECK) {
        fprintf(stderr, "Into verify %llu \n", NUM_INTO_VERIFY);
        fprintf(stderr, "Into 2nd check %llu, 1st check accepts %f%% pairs \n",
                NUM_INTO_2nd_CHECK, (double) NUM_INTO_2nd_CHECK / NUM_INTO_VERIFY * 100);
        fprintf(stderr, "Into extend %llu, 2nd check accepts %f%% pairs\n",
                NUM_INTO_EXTEND, (double) NUM_INTO_EXTEND / NUM_INTO_2nd_CHECK * 100);
        fprintf(stderr, "Into output %llu, extend accepts %f%% pairs\n\n",
                NUM_INTO_OUTPUT, (double) NUM_INTO_OUTPUT / NUM_INTO_EXTEND * 100);
        return 0;
    } else return 0;
}

int efficiency_check() {
    if(EFFICIENCY_CHECK) {
        fprintf(stderr, "Time before verify %.2fs, time in verify %.2fs \n"
        "  in verify takes %.2f%% of join time \n",
        TIME_BEFORE_VERIFY, TIME_IN_VERIFY, TIME_IN_VERIFY/(TIME_BEFORE_VERIFY+TIME_IN_VERIFY)*100);
        fprintf(stderr, "Time 1st check %.2fs \n"
        "  %.2f%% of verify time \n",
        TIME_1st_CHECK, TIME_1st_CHECK/TIME_IN_VERIFY*100);
        fprintf(stderr, "Time 2nd check %.2fs\n"
        "  %.2f%% of verify time \n",
        TIME_2nd_CHECK, TIME_2nd_CHECK/TIME_IN_VERIFY*100);
        fprintf(stderr, "Time extending %.2fs\n"
        "  %.2f%% of verify time \n",
        TIME_EXTEND, TIME_EXTEND/TIME_IN_VERIFY*100);
        fprintf(stderr, "Time cal match %.2fs\n"
        "  %.2f%% of 2nd check time \n",
        TIME_CAL_MATCH, TIME_CAL_MATCH/TIME_2nd_CHECK*100);
        return 0;
    } else return 0;
}
