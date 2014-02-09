#include "Read.h"

#include <sstream>
#include <string>

Read::Read(std::string line) {
    std::stringstream ss(line);
    std::string flag, chr, rnext, tmp;
    int map_quality, mate_pos, insertion_size;
    ss >> seq_name >> flag >> chr >> pos >> map_quality
       >> cigar_str >> rnext >> mate_pos >> insertion_size
       >> seq >> base_quality >> tmp >> md_str;
    
    // parse CIGAR string
    cigars = CIGAR::from_cigar_string(cigar_str);
    
    // parse MD tag
    mds = MD::from_md_string(md_str);
    
    // set left_most_pos
    if (cigars[0].type == CIGAR::SOFTCLIP) {
        left_most_pos = pos - cigars[0].length;
    } else {
        left_most_pos = pos;
    }
    
    // get variants from CIGAR & MD
    variants = VariantFromSAM::from_cigar_and_md(cigars, mds, left_most_pos, seq);
}
