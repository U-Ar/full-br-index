#ifndef INCLUDED_PERMUTED_LCP_HPP
#define INCLUDED_PERMUTED_LCP_HPP

#include "definitions.hpp"
#include "sparse_sd_vector.hpp"


namespace bri {

template<class sparse_bitvector_t = sparse_sd_vector>
class permuted_lcp 
{
public:
    permuted_lcp() {}

    /*
     * constructor with cache_config.
     * assume LCP and ISA are cached by the config.
     */
    permuted_lcp(sdsl::cache_config& cc)
    {
        sdsl::int_vector<> isa;
        sdsl::load_from_file(isa, sdsl::cache_file_name(sdsl::conf::KEY_ISA, cc));
        
        sdsl::int_vector<> lcp;
        sdsl::load_from_file(lcp, sdsl::cache_file_name(sdsl::conf::KEY_LCP, cc));
       
        n = isa.size();

        std::vector<bool> S(2*n+1,false);
        S[0] = true;

        ulint isa_val = isa[0];
        ulint lcp_prev = lcp[isa_val];
        ulint lcp_now;
        ulint pos = lcp_prev + 1;

        for (ulint i = 1; i < n; ++i)
        {
            S[pos] = true;
            isa_val = isa[i];
            lcp_now = lcp[isa_val];
            pos += lcp_now + 2;
            pos -= lcp_prev;
            lcp_prev = lcp_now;
        }
        S[pos] = true;

        u = pos + 1;
        std::vector<bool> ones_bv(u,false);
        std::vector<bool> zeros_bv(u,false);

        // generate run length encoding ones&zeros from S
        // S: unary encoding of j + PLCP[j]
        {
            bool bit_1 = true; // tell if we are counting 1 bits
            ulint cont_1 = 0, cont_0 = 0; // continuous number of 1 bits & 0 bits

            for (ulint i = 0; i < u; ++i)
            {
                if (S[i])
                {
                    if (bit_1)
                    {
                        cont_1++;
                    } else {
                        bit_1 = true;
                        cont_1++;
                        zeros_bv[cont_0-1] = true;
                    }
                } 
                else 
                {
                    if (bit_1)
                    {
                        bit_1 = false;
                        cont_0++;
                        ones_bv[cont_1-1] = true;
                    } else {
                        cont_0++;
                    }
                }
            }
            if (bit_1) ones_bv[cont_1-1] = true;
            else zeros_bv[cont_0-1] = true;
        }

        ones = sparse_bitvector_t(ones_bv);
        zeros = sparse_bitvector_t(zeros_bv);


    }

    /*
     * get PLCP[i]
     */
    ulint operator[](size_t i)
    {
        assert(i < n);
        ulint rank_0 = ones.rank(i+1);
        if (rank_0 > 0)
        {
            size_t res = zeros.select(rank_0-1) - i + 1;
            return res;
        }
        return 0;
    }

    ulint get_value(size_t i)
    {
        return operator[](i);
    }

    /*
     * serialize PLCP to the ostream
     */
    ulint serialize(std::ostream& out)
    {
        ulint w_bytes = 0;

        out.write((char*)&n,sizeof(n));
        out.write((char*)&u,sizeof(u));

        w_bytes += sizeof(n) + sizeof(u);

        if (n == 0) return w_bytes;

        w_bytes += ones.serialize(out);
        w_bytes += zeros.serialize(out);

        return w_bytes;
    }

    /*
     * load PLCP from the istream
     */
    void load(std::istream& in)
    {
        in.read((char*)&n, sizeof(n));
        in.read((char*)&u, sizeof(u));

        if (n == 0) return;

        ones.load(in);
        zeros.load(in);
    }

    ulint print_space()
    {
        ulint tot_bytes = 0;

        std::ofstream out("/dev/null");
        std::cout << "space for permuted_lcp:" << std::endl;

        auto bytes = ones.serialize(out);
        tot_bytes += bytes;
        std::cout << "- bitvector ones: " << bytes << " bytes" << std::endl;

        bytes = zeros.serialize(out);
        tot_bytes += bytes;
        std::cout << "- bitvector zeros: " << bytes << " bytes" << std::endl;

        return tot_bytes;
    }

    ulint get_space()
    {
        ulint tot_bytes = 0;

        std::ofstream out("/dev/null");

        tot_bytes += ones.serialize(out);
        tot_bytes += zeros.serialize(out);

        return tot_bytes;
    }

    ulint size()
    {
        return n;
    }

private:

    // length of LCP
    ulint n = 0;

    // length of ones & zeros
    ulint u = 0;

    // run length encoded PLCP representation
    sparse_bitvector_t ones;
    sparse_bitvector_t zeros;

};

};

#endif /* INCLUDED_PERMUTED_LCP_HPP */