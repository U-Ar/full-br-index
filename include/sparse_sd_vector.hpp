/*
 * sparse_sd_vector: a wrapper on sd_vector of the sdsl library, with support for rank/select1
 */

#ifndef INCLUDED_SPARSE_SD_VECTOR_HPP
#define INCLUDED_SPARSE_SD_VECTOR_HPP

#include "definitions.hpp"

namespace bri {

class sparse_sd_vector{
public:
    
    // empty contructor. all bits are 0
    sparse_sd_vector();

    // constructor. build using std::vector<bool>
    sparse_sd_vector(std::vector<bool> &b, bool enable_rank=true, bool enable_select=true);

    // constructor. build using bit_vector
    sparse_sd_vector(sdsl::bit_vector& bv, bool enable_rank=true, bool enable_select=true);

    // constructor. build from const iterators of positional vector
    template<class itr_t>
    sparse_sd_vector(const itr_t begin, const itr_t end) {
        sdv = sdsl::sd_vector<>(begin,end);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);

        u = sdv.size();
    }

    sparse_sd_vector& operator=(const sparse_sd_vector& other);

    inline bool rank_supported() { return rank_enabled; }
    inline bool select_supported() { return select_enabled; }

    /*
     * argument: position i 
     * returns: bit in position i
     * ACCESS ONLY
     */
    inline bool operator[](size_t i) { assert(i < size()); return sdv[i]; }

    inline bool at(size_t i) { return operator[](i); }

    /*
     * argument: position i 
     * returns: number of 1-bits in sdv[0...i-1]
     */
    inline ulint rank(size_t i)
    {
        assert(rank_enabled);
        assert(i <= size());
        return rank1(i);
    }

    /*
	 * argument: position 0<=i<=n
	 * returns: predecessor of i (position i excluded)
	 */
    inline size_t predecessor(size_t i)
    {
        assert(rank_enabled);
        assert(select_enabled);
        assert(rank(i)>0);
        return select(rank(i) - 1);
    }

    /*
	 * argument: position 0<=i<=n
	 * returns: rank of predecessor of i (i excluded) in
	 * bitvector space. If i does not have a predecessor,
	 * return rank of the last bit set in the bitvector
	 */
    inline ulint predecessor_rank_circular(ulint i)
    {
        return rank(i) == 0 ? number_of_1() - 1 : rank(i) - 1;
    }
    
    /*
	 * argument: bit number 0 <= i < m
	 * returns: gap between the i-th bit and (i-1)-th bit
	 */
    inline ulint gap_at(size_t i)
    {
        assert(select_enabled);
        assert(i<number_of_1());
        return i == 0 ? select(0)+1 : select(i)-select(i-1);
    }

    /*
	 * argument: ulint i >= 0
	 * returns: position of the i-th 1-bit
	 */
    inline size_t select(ulint i)
    {
        assert(select_enabled);
        assert(i<number_of_1());
        return select1(i+1);
    }

    /*
     * returns: size of the bitvector
     */
    inline ulint size() { return u; }

    /*
     * returns: number of 1s in the bitvector
     */
    inline ulint number_of_1() { 
        assert(rank_enabled);
        return rank1(u); 
    }

    /*
     * argument: ostream
     * returns: number of bytes written to ostream
     */
    ulint serialize(std::ostream& out);

    /*
     * load bitvector from istream
     * argument: istream
     */
    void load(std::istream& in);

private:

    //length of bitvector
    ulint u = 0;

    bool rank_enabled = true;
    bool select_enabled = true;

    sdsl::sd_vector<> sdv;
    sdsl::sd_vector<>::rank_1_type rank1;
    sdsl::sd_vector<>::select_1_type select1;

};

}; // namespace bri

#endif /* INCLUDED_SPARSE_SD_VECTOR_HPP */