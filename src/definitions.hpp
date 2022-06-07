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

};

#endif /* INCLUDED_DEFINITIONS_HPP */