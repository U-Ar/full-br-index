#include "huffman_string.hpp"

namespace bri {

/*
 * constructor.
 * takes a string as input and construct wt_huff on it
 */
huffman_string::huffman_string(std::string& s)
{
    s.push_back(0);
    sdsl::construct_im(wt, s.c_str(), 1);

    assert(wt.size()==s.size()-1);
}


}; // namespace bri