#include "iutest.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

#include "permuted_lcp.hpp"

using namespace bri;
using namespace sdsl;

IUTEST(PermutedLcpTest, Construction)
{
    std::string s("hello, world \n nice to meet you");

    cache_config cc;
    int_vector<8> text(s.size());
    for (ulint i = 0; i < s.size(); ++i)
            text[i] = (uchar)s[i];
    
    append_zero_symbol(text);

    store_to_cache(text, conf::KEY_TEXT, cc);
    construct_config::byte_algo_sa = SE_SAIS;

    construct_sa<8>(cc);
    construct_isa(cc);
    construct_lcp_kasai<8>(cc);

    permuted_lcp<> plcp(cc);

    remove(cache_file_name(conf::KEY_TEXT, cc));
    remove(cache_file_name(conf::KEY_SA, cc));
    remove(cache_file_name(conf::KEY_ISA, cc));
}

IUTEST(PermutedLcpTest, A10TextOperations)
{
    std::string s("aaaaaaaaaa");

    cache_config cc;
    int_vector<8> text(s.size());
    for (ulint i = 0; i < s.size(); ++i)
        text[i] = (uchar)s[i];
    append_zero_symbol(text);

    store_to_cache(text, conf::KEY_TEXT, cc);
    construct_config::byte_algo_sa = SE_SAIS;

    construct_sa<8>(cc);
    construct_isa(cc);
    construct_lcp_kasai<8>(cc);

    permuted_lcp<> plcp(cc);
    IUTEST_ASSERT_EQ(11,plcp.size());

    for (ulint i = 0; i < plcp.size(); ++i)
    {
        if (i == plcp.size()-1)
        {
            IUTEST_ASSERT_EQ(0,plcp[i]);
            IUTEST_ASSERT_EQ(0,plcp.get_value(i));
        } else {
            IUTEST_ASSERT_EQ(9-i,plcp[i]);
            IUTEST_ASSERT_EQ(9-i,plcp.get_value(i));
        }
    }

    std::ofstream ofs("test-tmp/permuted_lcp_test.tmp");
    auto w_bytes = plcp.serialize(ofs);
    IUTEST_ASSERT_EQ(plcp.print_space(),plcp.get_space());
    std::cout << "PLCP: " << w_bytes << " bytes" << std::endl;
    ofs.close();

    std::ifstream ifs("test-tmp/permuted_lcp_test.tmp");
    plcp.load(ifs);
    ifs.close();

    IUTEST_ASSERT_EQ(11,plcp.size());
    for (ulint i = 0; i < plcp.size(); ++i)
    {
        if (i == plcp.size()-1)
        {
            IUTEST_ASSERT_EQ(0,plcp[i]);
            IUTEST_ASSERT_EQ(0,plcp.get_value(i));
        } else {
            IUTEST_ASSERT_EQ(9-i,plcp[i]);
            IUTEST_ASSERT_EQ(9-i,plcp.get_value(i));
        }
    }

    remove(cache_file_name(conf::KEY_TEXT, cc));
    remove(cache_file_name(conf::KEY_SA, cc));
    remove(cache_file_name(conf::KEY_ISA, cc));
}

IUTEST(PermutedLcpTest, TwoCharTextOperations)
{
    std::string s("aaaaabbbbbaaaaa");

    cache_config cc;
    int_vector<8> text(s.size());
    for (ulint i = 0; i < s.size(); ++i)
        text[i] = (uchar)s[i];
    append_zero_symbol(text);

    store_to_cache(text, conf::KEY_TEXT, cc);
    construct_config::byte_algo_sa = SE_SAIS;

    construct_sa<8>(cc);
    construct_isa(cc);
    construct_lcp_kasai<8>(cc);

    permuted_lcp<> plcp(cc);
    IUTEST_ASSERT_EQ(16,plcp.size());
    
    IUTEST_ASSERT_EQ(5,plcp[0]);
    IUTEST_ASSERT_EQ(4,plcp[1]);
    IUTEST_ASSERT_EQ(3,plcp[2]);
    IUTEST_ASSERT_EQ(2,plcp[3]);
    IUTEST_ASSERT_EQ(1,plcp[4]);
    IUTEST_ASSERT_EQ(4,plcp[5]);
    IUTEST_ASSERT_EQ(3,plcp[6]);
    IUTEST_ASSERT_EQ(2,plcp[7]);
    IUTEST_ASSERT_EQ(1,plcp[8]);
    IUTEST_ASSERT_EQ(0,plcp[9]);
    IUTEST_ASSERT_EQ(4,plcp[10]);
    IUTEST_ASSERT_EQ(3,plcp[11]);
    IUTEST_ASSERT_EQ(2,plcp[12]);
    IUTEST_ASSERT_EQ(1,plcp[13]);
    IUTEST_ASSERT_EQ(0,plcp[14]);
    IUTEST_ASSERT_EQ(0,plcp[15]);

    std::ofstream ofs("test-tmp/permuted_lcp_test2.tmp");
    auto w_bytes = plcp.serialize(ofs);
    IUTEST_ASSERT_EQ(plcp.print_space(),plcp.get_space());
    std::cout << "PLCP: " << w_bytes << " bytes" << std::endl;
    ofs.close();

    std::ifstream ifs("test-tmp/permuted_lcp_test2.tmp");
    plcp.load(ifs);
    ifs.close();

    IUTEST_ASSERT_EQ(16,plcp.size());
    IUTEST_ASSERT_EQ(5,plcp[0]);
    IUTEST_ASSERT_EQ(4,plcp[1]);
    IUTEST_ASSERT_EQ(3,plcp[2]);
    IUTEST_ASSERT_EQ(2,plcp[3]);
    IUTEST_ASSERT_EQ(1,plcp[4]);
    IUTEST_ASSERT_EQ(4,plcp[5]);
    IUTEST_ASSERT_EQ(3,plcp[6]);
    IUTEST_ASSERT_EQ(2,plcp[7]);
    IUTEST_ASSERT_EQ(1,plcp[8]);
    IUTEST_ASSERT_EQ(0,plcp[9]);
    IUTEST_ASSERT_EQ(4,plcp[10]);
    IUTEST_ASSERT_EQ(3,plcp[11]);
    IUTEST_ASSERT_EQ(2,plcp[12]);
    IUTEST_ASSERT_EQ(1,plcp[13]);
    IUTEST_ASSERT_EQ(0,plcp[14]);
    IUTEST_ASSERT_EQ(0,plcp[15]);

    remove(cache_file_name(conf::KEY_TEXT, cc));
    remove(cache_file_name(conf::KEY_SA, cc));
    remove(cache_file_name(conf::KEY_ISA, cc));
}


IUTEST(PermutedLcpTest, A1000000TextOperations)
{
    std::string s;
    for (ulint i = 0; i < 1000000; ++i) s.push_back('a');

    cache_config cc;
    int_vector<8> text(s.size());
    for (ulint i = 0; i < s.size(); ++i)
        text[i] = (uchar)s[i];
    append_zero_symbol(text);

    store_to_cache(text, conf::KEY_TEXT, cc);
    construct_config::byte_algo_sa = SE_SAIS;

    construct_sa<8>(cc);
    construct_isa(cc);
    construct_lcp_kasai<8>(cc);

    permuted_lcp<> plcp(cc);
    IUTEST_ASSERT_EQ(1000001,plcp.size());

    for (ulint i = 0; i < plcp.size(); ++i)
    {
        if (i == plcp.size()-1)
        {
            IUTEST_ASSERT_EQ(0,plcp[i]);
            IUTEST_ASSERT_EQ(0,plcp.get_value(i));
        } else {
            IUTEST_ASSERT_EQ(999999-i,plcp[i]);
            IUTEST_ASSERT_EQ(999999-i,plcp.get_value(i));
        }
    }

    std::ofstream ofs("test-tmp/permuted_lcp_test.tmp");
    auto w_bytes = plcp.serialize(ofs);
    IUTEST_ASSERT_EQ(plcp.print_space(),plcp.get_space());
    std::cout << "PLCP: " << w_bytes << " bytes" << std::endl;
    ofs.close();

    std::ifstream ifs("test-tmp/permuted_lcp_test.tmp");
    plcp.load(ifs);
    ifs.close();

    IUTEST_ASSERT_EQ(1000001,plcp.size());
    for (ulint i = 0; i < plcp.size(); ++i)
    {
        if (i == plcp.size()-1)
        {
            IUTEST_ASSERT_EQ(0,plcp[i]);
            IUTEST_ASSERT_EQ(0,plcp.get_value(i));
        } else {
            IUTEST_ASSERT_EQ(999999-i,plcp[i]);
            IUTEST_ASSERT_EQ(999999-i,plcp.get_value(i));
        }
    }

    remove(cache_file_name(conf::KEY_TEXT, cc));
    remove(cache_file_name(conf::KEY_SA, cc));
    remove(cache_file_name(conf::KEY_ISA, cc));
}



/* test using external textfile
IUTEST(PermutedLcpTest, BigText)
{
    std::ifstream ifs("dna.50MB");
    std::stringstream buf;
    buf << ifs.rdbuf();
    std::string s(buf.str());
    ifs >> s;
    std::cout << "testing DNA " << std::endl;
    std::cout << s.size() << std::endl;

    cache_config cc;
    int_vector<8> text(s.size());
    for (ulint i = 0; i < s.size(); ++i)
        text[i] = (uchar)s[i];
    append_zero_symbol(text);

    store_to_cache(text, conf::KEY_TEXT, cc);
    construct_config::byte_algo_sa = SE_SAIS;

    construct_sa<8>(cc);
    construct_isa(cc);

    permuted_lcp<> plcp(cc);

    construct_lcp_kasai<8>(cc);
    int_vector<> lcp;
    load_from_file(lcp,cache_file_name(conf::KEY_LCP,cc));
    int_vector<> isa;
    load_from_file(isa,cache_file_name(conf::KEY_ISA,cc));

    IUTEST_EXPECT_EQ(isa.size(),plcp.size());
    for (ulint i = 0; i < isa.size(); ++i)
    {
        IUTEST_EXPECT_EQ(lcp[isa[i]],plcp[i]);
    }


    remove(cache_file_name(conf::KEY_TEXT, cc));
    remove(cache_file_name(conf::KEY_SA, cc));
    remove(cache_file_name(conf::KEY_ISA, cc));
    remove(cache_file_name(conf::KEY_LCP, cc));
}
*/