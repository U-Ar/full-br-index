#include "sparse_sd_vector.hpp"

namespace bri {

// empty contructor. all bits are 0
sparse_sd_vector::sparse_sd_vector() : rank_enabled(true), select_enabled(true) {}

// constructor. build using std::vector<bool>
sparse_sd_vector::sparse_sd_vector(std::vector<bool> &b, bool enable_rank, bool enable_select) 
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

// constructor. build using bit_vector
sparse_sd_vector::sparse_sd_vector(sdsl::bit_vector& bv, bool enable_rank, bool enable_select)
{
    rank_enabled = enable_rank;
    select_enabled = enable_select;

    u = bv.size();

    sdv = sdsl::sd_vector<>(bv);
    if (rank_enabled) rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
    if (select_enabled) select1 = sdsl::sd_vector<>::select_1_type(&sdv);
}

sparse_sd_vector& sparse_sd_vector::operator=(const sparse_sd_vector& other)
{
    u = other.sdv.size();
    if (other.rank_enabled) rank_enabled = true;
    if (other.select_enabled) select_enabled = true;

    sdv = sdsl::sd_vector<>(other.sdv);
    if (rank_enabled) rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
    if (select_enabled) select1 = sdsl::sd_vector<>::select_1_type(&sdv);

    return *this;
}

/*
 * serialize the bitvector to the ostream
 * argument: ostream
 * returns: number of bytes written to ostream
 */
ulint sparse_sd_vector::serialize(std::ostream& out) 
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
void sparse_sd_vector::load(std::istream& in) 
{
    in.read((char*)&u, sizeof(u));

    if (u == 0) return;

    sdv.load(in);
    if (rank_enabled) rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
    if (select_enabled) select1 = sdsl::sd_vector<>::select_1_type(&sdv);
}


}; // namespace bri