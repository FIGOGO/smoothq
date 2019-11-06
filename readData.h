#pragma once

#include "util.h"
#include <queue>
#include <fstream>

inline std::string revstring(const std::string& input)
{
	std::string output;
	output.resize(input.size());

	transform(input.rbegin(), input.rend(), output.begin(),
	        [](char c) -> char {return trans[c];}
	        );

	return output;
};


inline void add_read(const std::string& read, const std::string& tag, const int index) {
	if (read.size() <= MIN_LEN) return;

    // add the forward sequence
    Data newdata;
    newdata.str = read;
    newdata.index = index;
    newdata.tag = tag;
    newdata.len = read.size();
    newdata.fwrv = 0;
    newdata.cstr = newdata.str.c_str();
    oridata.push_back(newdata);

    // add the reverse sequence
    newdata.str = revstring(read);
    newdata.fwrv = 1;
    newdata.cstr = newdata.str.c_str();
    oridata.push_back(newdata);
};


//extract reads from fasta file
void readfasta(const std::string& filename) {
    std::ifstream input_stream(filename);

    std::string line;
    std::string read = "\n";
    std::queue<std::string> tag;
    int index = 1;

    while (getline(input_stream, line)) {
        if (line.empty())
            continue;
        else if (line[0] == '>') {
            tag.push(line.substr(1, line.find(" ")));
            read.erase(remove(read.begin(), read.end(), '\n'), read.end());
            if (!read.empty()) {
                add_read(read, tag.front(), index);
                index++;
                tag.pop();
                read.clear();
            }
            else {
                read = "\n";
            }
        }
        else {
            read += line;
        }
    }
    // add last read
    if (!read.empty()) {
        add_read(read, tag.front(), index);
		tag.pop();
    }
    for (int i = 0; i < oridata.size(); i++) {
        oridata[i].cstr = oridata[i].str.c_str();
    }

    fprintf(stderr, "Load %d sequences\n", oridata.size());
};

