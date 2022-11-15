#ifndef INCLUDED_DEFINITIONS_HPP
#define INCLUDED_DEFINITIONS_HPP

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <sdsl/construct.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/lcp_bitcompressed.hpp>

namespace bri {

typedef uint64_t ulint;
typedef long int lint;
typedef uint32_t uint;
typedef unsigned char uchar;

typedef std::pair<ulint, ulint> range_t;

struct range_hash {
    std::size_t operator() (range_t const& range) const;
};

// sample maintained during the search
struct br_sample {
    /*
     * state variables for left_extension & right_extension
     * range: SA range of P
     * j: SA[p]
     * d: offset between starting position of the pattern & j
     * rangeR: correspondents to range in reversed text
     * len: current pattern length
     */
    range_t range, rangeR;
    ulint j, d, len;

    br_sample() {}

    br_sample(range_t range_, 
              range_t rangeR_,
              ulint j_,
              ulint d_,
              ulint len_) 
              :
              range(range_),
              rangeR(rangeR_),
              j(j_),
              d(d_),
              len(len_) {}

    inline void set_values(range_t range, 
                    range_t rangeR,
                    ulint j,
                    ulint d,
                    ulint len)
    {
        this->range=range; this->rangeR=rangeR;
        this->j=j;this->d=d;this->len=len;
    }
    
    // the pattern does not exist
    inline bool is_invalid() const
    {
        return (range.first > range.second) || (rangeR.first > rangeR.second);
    }

    // the pattern exists
    inline bool is_valid() const 
    {
        return (range.first <= range.second) && (rangeR.first <= rangeR.second);
    }

    inline bool is_leaf() const
    {
        return size() == 1;
    }

    inline bool contains(br_sample const& other) const
    {
        return range.first <= other.range.first && other.range.second <= range.second;
    }

    // range size
    inline ulint size() const
    {
        return is_valid() ? (range.second + 1 - range.first) : 0;
    }
    
};

std::string get_time(ulint time);
uchar bitsize(ulint x);
void header_error();
ulint get_number_of_patterns(std::string header);
ulint get_patterns_length(std::string header);

};

#endif /* INCLUDED_DEFINITIONS_HPP */