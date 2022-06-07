/*
 * huffman_string: a wrapper on wt_huff of the sdsl library, with support for rank/select
 */

#ifndef INCLUDED_HUFFMAN_STRING_HPP
#define INCLUDED_HUFFMAN_STRING_HPP

#include "definitions.hpp"

namespace bri {

class huffman_string {

public:
    huffman_string() {}

    /*
     * constructor.
     * takes a string as input and construct wt_huff on it
     */
    huffman_string(std::string& s)
    {
        s.push_back(0);
        sdsl::construct_im(wt, s.c_str(), 1);

        assert(wt.size()==s.size()-1);
    }

    uchar operator[](size_t i) 
    {

        assert(i<wt.size());
        return wt[i];

    }

    /*
     * length of the string
     */
    size_t size()
    {
        return wt.size();
    }

    /*
     * number of character c in S[0...i-1]
     */
    ulint rank(size_t i, uchar c)
    {
        assert(i<=wt.size());
        return wt.rank(i,c);
    }

    /*
     * position of i-th character c
     * i starts from 0
     */
    ulint select(ulint i, uchar c)
    {
        return wt.select(i+1,c);
    }

    /*
     * serialize the index to the ostream
     */
    ulint serialize(std::ostream& out)
    {
        return wt.serialize(out);
    }

    /*
     * load the index from the ostream
     */
    void load(std::istream& in)
    {
        wt.load(in);
    }

private:

    sdsl::wt_huff<> wt;
    
};

};

#endif /* INCLUDED_HUFFMAN_STRING_HPP */