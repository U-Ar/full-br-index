/*
 * sparse_sd_vector: a wrapper on sd_vector of the sdsl library, with support for rank/select1
 */

#ifndef INCLUDED_SPARSE_SD_VECTOR_HPP
#define INCLUDED_SPARSE_SD_VECTOR_HPP

#include "definitions.hpp"

#ifndef ulint
typedef uint64_t ulint;
#endif

#ifndef uint
typedef uint32_t uint;
#endif 

namespace bri {

class sparse_sd_vector{
public:

    /*
     * empty contructor. all bits are 0
     */
    sparse_sd_vector(bool enable_rank=true, bool enable_select=true) :
        rank_enabled(enable_rank), select_enabled(enable_select)
    {}

    /*
     * constructor. build using std::vector<bool>
     */
    sparse_sd_vector(std::vector<bool> &b, bool enable_rank=true, bool enable_select=true)
    {
        if (b.size() == 0) return;

        rank_enabled = enable_rank;
        select_enabled = enable_select;

        u = b.size();

        sdsl::bit_vector bv(u);
        for (size_t i = 0; i < u; ++i) bv[i] = b[i];

        sdv = sdsl::sd_vector<>(bv);
        if (rank_enabled) rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        if (select_enabled) select1 = sdsl::sd_vector<>::select_1_type(&sdv);

    }

    /*
     * constructor. build using bit_vector
     */
    sparse_sd_vector(sdsl::bit_vector& bv, bool enable_rank=true, bool enable_select=true)
    {
        rank_enabled = enable_rank;
        select_enabled = enable_select;

        u = bv.size();

        sdv = sdsl::sd_vector<>(bv);
        if (rank_enabled) rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        if (select_enabled) select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    /*
     * substitution operator.
     */
    sparse_sd_vector& operator=(const sparse_sd_vector& other)
    {
        u = other.sdv.size();
        if (other.rank_enabled) rank_enabled = true;
        if (other.select_enabled) select_enabled = true;

        sdv = sdsl::sd_vector<>(other.sdv);
        if (rank_enabled) rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        if (select_enabled) select1 = sdsl::sd_vector<>::select_1_type(&sdv);

        return *this;
    }

    bool rank_supported() { return rank_enabled; }
    bool select_supported() { return select_enabled; }

    /*
     * argument: position i 
     * returns: bit in position i
     * ACCESS ONLY
     */
    bool operator[](size_t i) 
    {
        assert(i < size());
        return sdv[i];
    }

    bool at(size_t i)
    {
        return operator[](i);
    }

    /*
     * argument: position i 
     * returns: number of 1-bits in sdv[0...i-1]
     */
    ulint rank(size_t i)
    {
        assert(rank_enabled);
        assert(i <= size());

        return rank1(i);
    }

    /*
	 * argument: position 0<=i<=n
	 * returns: predecessor of i (position i excluded)
	 */
    size_t predecessor(size_t i)
    {
        assert(rank_enabled);
        assert(select_enabled);
        assert(rank(i)>0);

        return select(rank(i) - 1);
    }

    /*
	 * input: position 0<=i<=n
	 * output: rank of predecessor of i (i excluded) in
	 * bitvector space. If i does not have a predecessor,
	 * return rank of the last bit set in the bitvector
	 */
    ulint predecessor_rank_circular(ulint i)
    {
        return rank(i) == 0 ? number_of_1() - 1 : rank(i) - 1;
    }
    
    ulint gap_at(size_t i)
    {
        assert(select_enabled);
        assert(i<number_of_1());
        
        if (i == 0) return select(0) + 1;

        return select(i) - select(i-1);
    }

    /*
	 * argument: ulint i >= 0
	 * returns: position of the i-th 1-bit
	 */
    size_t select(ulint i)
    {
        assert(select_enabled);
        assert(i<number_of_1());
        
        return select1(i+1);
    }

    /*
     * returns: size of the bitvector
     */
    ulint size() { return u; }

    /*
     * returns: number of 1s in the bitvector
     */
    ulint number_of_1() { 

        assert(rank_enabled);
        return rank1(u); 

    }

    /*
     * argument: ostream
     * returns: number of bytes written to ostream
     */
    ulint serialize(std::ostream& out)
    {
        ulint w_bytes = 0;

        out.write((char*)&u, sizeof(u));
        w_bytes += sizeof(u);

        if (u == 0) return w_bytes;

        w_bytes += sdv.serialize(out);
        return w_bytes;
    }

    /*
     * load bitvector from istream
     * argument: istream
     */
    void load(std::istream& in)
    {
        in.read((char*)&u, sizeof(u));

        if (u == 0) return;

        sdv.load(in);
        if (rank_enabled) rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        if (select_enabled) select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

private:

    //length of bitvector
    ulint u = 0;

    bool rank_enabled = true;
    bool select_enabled = true;

    sdsl::sd_vector<> sdv;
    sdsl::sd_vector<>::rank_1_type rank1;
    sdsl::sd_vector<>::select_1_type select1;

};

};

#endif /* INCLUDED_SPARSE_SD_VECTOR_HPP */