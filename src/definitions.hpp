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

struct range_hash
{
    std::size_t operator() (range_t const& range) const
    {
        auto hash1 = std::hash<ulint>()(range.first);
        auto hash2 = std::hash<ulint>()(range.second);
        std::size_t seed = 0;
        seed ^= hash1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }

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
    
    br_sample(): range(), rangeR() {}

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

    void set_values(range_t range_, 
                    range_t rangeR_,
                    ulint j_,
                    ulint d_,
                    ulint len_)
    {
        range = range_;
        rangeR = rangeR_;
        j = j_;
        d = d_;
        len = len_;
    }

    // the pattern does not exist
    inline bool is_invalid() const
    {
        return (range.first > range.second) || (rangeR.first > rangeR.second);
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
        return range.second + 1 - range.first;
    }
};

};

#endif /* INCLUDED_DEFINITIONS_HPP */