#include "smoothq_implementations.h"
using namespace std;

int smoothq_vec() {
	double time_read_cpu, time_generate_signature_cpu, time_generate_table_cpu,
			time_join2_cpu;
	double time_read, time_generate_signature, time_generate_table, time_join2;
	double current_time_cpu, current_time;

	current_time = time(NULL);
	current_time_cpu = get_cpu_time();
	// Generate signatures
	srand((unsigned int) current_time);
	vector<Signature_inter> ht = generate_signature_vec();
	time_generate_signature_cpu=get_cpu_time()-current_time_cpu;
	time_generate_signature=time(NULL)-current_time;

	current_time = time(NULL);
	current_time_cpu = get_cpu_time();
	// Generate hash table for candidates
	vector<Signature_inter> htable = generate_table_vec(ht);
	time_generate_table_cpu=get_cpu_time()-current_time_cpu;
    time_generate_table=time(NULL)-current_time;

	current_time = time(NULL);
	current_time_cpu = get_cpu_time();
	// Generate hash table for candidates
	join(htable);
	cout << flush;
	fprintf(stderr, "Join takes %f s cpu \n",
			(time_join2_cpu=get_cpu_time()-current_time_cpu));
	fprintf(stderr, "     takes %f s real clock \n",
			(time_join2=time(NULL)-current_time));


	fprintf(stderr, "Main algorithm total time cpu: %f s, real clock %f s \n",
			time_read_cpu+time_generate_signature_cpu+time_generate_table_cpu+time_join2_cpu,
			time_read+time_generate_signature+time_generate_table+time_join2);
	return 0;
}

int parameter_parsing(int argc, char **argv) {
    const char *opt_str = "f:r:t:c:q:e:m:o:i:";
    int opt;

    while ((opt = getopt (argc, argv, opt_str)) != -1)
        switch (opt)
        {
            case 'f':
                TOPK = stof(optarg);
                break;
            case 'c':
                THRESHOLD = atoi(optarg);
                break;
            case 'r':
                SAMPLING_RATE = stof(optarg);
                break;
            case 't':
                NUM_THREADS = atoi(optarg);
                break;
            case 'q':
                SIZE_Q = atoi(optarg);
                break;
            case 'e':
                SIZE_CGK = atoi(optarg);
                break;
            case 'm':
                SIZE_SMQ = atoi(optarg);
                break;
            case 'o':
                OUT_FORMAT = optarg;
                cout << OUT_FORMAT << endl;
                if ((OUT_FORMAT!="m4") && (OUT_FORMAT!="paf")) {
                    throw std::invalid_argument("Please specify a valid output format: m4 or paf");
                }
                break;
            case 'i':
                OPT_INPUT = optarg;
                break;
            case '?':
                if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                             "Unknown option character `\\x%x'.\n",
                             optopt);
                return 1;
            default:
                continue;
        }
    return 0;
}

int self_join(string filename)
{
	double starting_time_cpu = get_cpu_time();
	double starting_time = time(NULL);

	// Read data from input file
	readfasta(filename);
	NUM_READS_DATA1 = oridata.size();
	fprintf(stderr, "Finished reading input file in %f s cpu \n",
			(get_cpu_time()-starting_time_cpu));
	fprintf(stderr, "                            in %f s real clock \n",
			(time(NULL)-starting_time));

	// Print out parameters used in the program
    print_working_env();

    // smoothq algorithm with vector based hash table
	smoothq_vec();

	fprintf(stderr, "Total running time cpu: %f s, real clock %f s \n",
	        get_cpu_time()-starting_time_cpu, time(NULL)-starting_time);

	// By default the check is off
    number_check();
    efficiency_check();

	return 0;
}

int two_join(std::string filename) {
    double starting_time_cpu = get_cpu_time();
    double starting_time = time(NULL);

    // Read data from input file
    readfasta(filename);
    NUM_READS_DATA1 = oridata.size();
    readfasta(OPT_INPUT);
    NUM_READS_DATA2 = oridata.size() - NUM_READS_DATA1;
    if (NUM_READS_DATA1 == 0 || NUM_READS_DATA2 == 0) {
        throw std::range_error("Input dataset is empty");
    }
    SINGULAR = 2;

    fprintf(stderr, "total size: %d; num reads 1: %d; num reads 2: %d\n",
            oridata.size(), NUM_READS_DATA1, NUM_READS_DATA2);
    fprintf(stderr, "Finished reading input file in %f s cpu \n",
            (get_cpu_time()-starting_time_cpu));
    fprintf(stderr, "                            in %f s real clock \n",
            (time(NULL)-starting_time));

    // Print out parameters used in the program
    print_working_env();

    // smoothq algorithm with vector based hash table
    smoothq_vec();

    fprintf(stderr, "Total running time cpu: %f s, real clock %f s \n",
            get_cpu_time()-starting_time_cpu, time(NULL)-starting_time);

    return 0;
}

int main(int argc, char **argv) {
    if (argc == 0) {
        fprintf(stderr, "ERROR: Please specify the input file.\n");
        return -1;
    }

    string filename = argv[1];
    parameter_parsing(argc, argv);

    if (OPT_INPUT == "") {
        fprintf(stderr, "Working on data: %s \n", filename.c_str());
        self_join(filename);
    } else {
        fprintf(stderr, "Working on datasets: %s and %s \n",
                filename.c_str(), OPT_INPUT.c_str());
        two_join(filename);
    }
    return 0;
}
