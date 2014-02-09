#ifndef __Read__
#define __Read__

#include <iostream>
#include <string>
#include <vector>
#include "MD.h"
#include "CIGAR.h"
#include "VariantFromSAM.h"

class Read {
public:
    int pos, left_most_pos;
    std::string seq_name;
    std::string seq, base_quality, md_str, cigar_str;
    std::vector<MD> mds;
    std::vector<CIGAR> cigars;
    std::vector<VariantFromSAM> variants;
    Read(std::string line);

    friend std::ostream & operator<<(std::ostream &stream, const Read &read) {
        stream << read.seq_name << std::endl;
        stream << read.cigar_str << std::endl;
        for (int i = 0;i < read.cigars.size();i++) {
            stream << read.cigars[i] << ", ";
        }
        stream << std::endl << read.md_str << std::endl;
        for (int i = 0;i < read.mds.size();i++) {
            stream << read.mds[i] << ", ";
        }
        stream << std::endl << "variants:" << std::endl;
        for (int i = 0;i < read.variants.size();i++) {
            stream << read.variants[i] << ", ";
        }
        stream << std::endl;
        return stream;
    }
};

#endif
