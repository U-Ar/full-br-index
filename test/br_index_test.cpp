#include "iutest.hpp"
#include <vector>
#include <fstream>
#include <string>

#include "../src/br_index.hpp"
#include "../src/br_index_naive.hpp"

using namespace bri;

template<class T>
void print_vec(std::vector<T>& vec)
{
    for (size_t i = 0; i < vec.size(); ++i)
        std::cout << vec[i] << " ";
    std::cout << std::endl;
}

void print_range(range_t rg)
{
    std::cout << "1st val: " << rg.first <<
    " 2nd val: " << rg.second << std::endl;
}

IUTEST(BrIndexTest, BasicLocate)
{
    std::string s("aaaaaaaaaaaaaaaaaaaa");
    br_index<> idx(s);
    IUTEST_ASSERT_EQ(20,idx.text_size());
    IUTEST_ASSERT_EQ(21,idx.bwt_size());
    range_t range = idx.full_range();
    IUTEST_ASSERT_EQ(20,range.second);
    IUTEST_ASSERT_EQ(0,range.first);

    std::vector<ulint> res = idx.locate("a");
    IUTEST_ASSERT_EQ(20,res.size());
    std::sort(res.begin(),res.end());
    for (ulint i = 0; i < res.size(); ++i)
    {
        IUTEST_ASSERT_EQ(i,res[i]);
    }

    res = idx.locate("aa");
    IUTEST_ASSERT_EQ(19,res.size());
    std::sort(res.begin(),res.end());
    for (ulint i = 0; i < res.size(); ++i)
    {
        IUTEST_ASSERT_EQ(i,res[i]);
    }

    res = idx.locate("aa",false);
    IUTEST_ASSERT_EQ(19,res.size());
    std::sort(res.begin(),res.end());
    for (ulint i = 0; i < res.size(); ++i)
    {
        IUTEST_ASSERT_EQ(i,res[i]);
    }
}

IUTEST(BrIndexTest, PeriodicTextLocate)
{
    std::string s("abcdabcdabcdabcdhello");
    br_index<> idx(s);
    auto vec = idx.locate("abcd");
    std::sort(vec.begin(),vec.end());
    IUTEST_ASSERT_EQ(0,vec[0]);
    IUTEST_ASSERT_EQ(4,vec[1]);
    IUTEST_ASSERT_EQ(8,vec[2]);
    IUTEST_ASSERT_EQ(12,vec[3]);
    vec = idx.locate("abcd",false);
    std::sort(vec.begin(),vec.end());
    IUTEST_ASSERT_EQ(0,vec[0]);
    IUTEST_ASSERT_EQ(4,vec[1]);
    IUTEST_ASSERT_EQ(8,vec[2]);
    IUTEST_ASSERT_EQ(12,vec[3]);

}

IUTEST(BrIndexTest, PhiPhiI)
{
    std::string s("abcdabcdabcdabcdhello");
    br_index<> idx(s);
    IUTEST_EXPECT_EQ(21,idx.Phi(0));
    IUTEST_EXPECT_EQ(0,idx.Phi(4));
    IUTEST_EXPECT_EQ(4,idx.Phi(8));
    IUTEST_EXPECT_EQ(8,idx.Phi(12));
    IUTEST_EXPECT_EQ(12,idx.Phi(1));
    IUTEST_EXPECT_EQ(1,idx.Phi(5));
    IUTEST_EXPECT_EQ(5,idx.Phi(9));
    IUTEST_EXPECT_EQ(9,idx.Phi(13));
    IUTEST_EXPECT_EQ(13,idx.Phi(2));
    IUTEST_EXPECT_EQ(2,idx.Phi(6));
    IUTEST_EXPECT_EQ(6,idx.Phi(10));
    IUTEST_EXPECT_EQ(10,idx.Phi(14));
    IUTEST_EXPECT_EQ(14,idx.Phi(3));
    IUTEST_EXPECT_EQ(3,idx.Phi(7));

    IUTEST_EXPECT_EQ(0,idx.PhiI(21));
    IUTEST_EXPECT_EQ(4,idx.PhiI(0));
    IUTEST_EXPECT_EQ(8,idx.PhiI(4));
    IUTEST_EXPECT_EQ(12,idx.PhiI(8));
    IUTEST_EXPECT_EQ(1,idx.PhiI(12));
    IUTEST_EXPECT_EQ(5,idx.PhiI(1));
    IUTEST_EXPECT_EQ(9,idx.PhiI(5));
    IUTEST_EXPECT_EQ(13,idx.PhiI(9));
    IUTEST_EXPECT_EQ(2,idx.PhiI(13));
    IUTEST_EXPECT_EQ(6,idx.PhiI(2));
    IUTEST_EXPECT_EQ(10,idx.PhiI(6));
    IUTEST_EXPECT_EQ(14,idx.PhiI(10));
    IUTEST_EXPECT_EQ(3,idx.PhiI(14));
    IUTEST_EXPECT_EQ(7,idx.PhiI(3));

}



IUTEST(BrIndexTest, DNALikeTextLocate)
{
    std::string s("AAAATGCCGCCGCCATAAA");
    br_index<> idx(s);

    auto vec = idx.locate("C");
    std::sort(vec.begin(),vec.end());
    print_vec<>(vec);
    IUTEST_EXPECT_EQ(6,vec[0]);
    IUTEST_EXPECT_EQ(7,vec[1]);
    IUTEST_EXPECT_EQ(9,vec[2]);
    IUTEST_EXPECT_EQ(10,vec[3]);
    IUTEST_EXPECT_EQ(12,vec[4]);
    IUTEST_EXPECT_EQ(13,vec[5]);

    vec = idx.locate("CC");
    std::sort(vec.begin(),vec.end());
    print_vec<>(vec);
    IUTEST_EXPECT_EQ(6,vec[0]);
    IUTEST_EXPECT_EQ(9,vec[1]);
    IUTEST_EXPECT_EQ(12,vec[2]);

    vec = idx.locate("CC",false);
    std::sort(vec.begin(),vec.end());
    print_vec<>(vec);
    IUTEST_EXPECT_EQ(6,vec[0]);
    IUTEST_EXPECT_EQ(9,vec[1]);
    IUTEST_EXPECT_EQ(12,vec[2]);


    vec = idx.locate("GCC");
    std::sort(vec.begin(),vec.end());
    print_vec<>(vec);
    IUTEST_EXPECT_EQ(5,vec[0]);
    IUTEST_EXPECT_EQ(8,vec[1]);
    IUTEST_EXPECT_EQ(11,vec[2]);

    vec = idx.locate("GCC",false);
    std::sort(vec.begin(),vec.end());
    print_vec<>(vec);
    IUTEST_EXPECT_EQ(5,vec[0]);
    IUTEST_EXPECT_EQ(8,vec[1]);
    IUTEST_EXPECT_EQ(11,vec[2]);
}



/*
 * below: tests for br_index_naive (same as tests for br_index)
 */


IUTEST(BrIndexNaiveTest, BasicLocate)
{
    std::string s("aaaaaaaaaaaaaaaaaaaa");
    br_index_naive<> idx(s);
    IUTEST_ASSERT_EQ(20,idx.text_size());
    IUTEST_ASSERT_EQ(21,idx.bwt_size());
    range_t range = idx.full_range();
    IUTEST_ASSERT_EQ(20,range.second);
    IUTEST_ASSERT_EQ(0,range.first);

    std::vector<ulint> res = idx.locate("a");
    IUTEST_ASSERT_EQ(20,res.size());
    std::sort(res.begin(),res.end());
    for (ulint i = 0; i < res.size(); ++i)
    {
        IUTEST_ASSERT_EQ(i,res[i]);
    }

    res = idx.locate("aa");
    IUTEST_ASSERT_EQ(19,res.size());
    std::sort(res.begin(),res.end());
    for (ulint i = 0; i < res.size(); ++i)
    {
        IUTEST_ASSERT_EQ(i,res[i]);
    }

    res = idx.locate("aa",false);
    IUTEST_ASSERT_EQ(19,res.size());
    std::sort(res.begin(),res.end());
    for (ulint i = 0; i < res.size(); ++i)
    {
        IUTEST_ASSERT_EQ(i,res[i]);
    }
}


IUTEST(BrIndexNaiveTest, PeriodicTextLocate)
{
    std::string s("abcdabcdabcdabcdhello");
    br_index_naive<> idx(s);
    auto vec = idx.locate("abcd");
    std::sort(vec.begin(),vec.end());
    IUTEST_ASSERT_EQ(0,vec[0]);
    IUTEST_ASSERT_EQ(4,vec[1]);
    IUTEST_ASSERT_EQ(8,vec[2]);
    IUTEST_ASSERT_EQ(12,vec[3]);

    vec = idx.locate("abcd",false);
    std::sort(vec.begin(),vec.end());
    IUTEST_ASSERT_EQ(0,vec[0]);
    IUTEST_ASSERT_EQ(4,vec[1]);
    IUTEST_ASSERT_EQ(8,vec[2]);
    IUTEST_ASSERT_EQ(12,vec[3]);

}

IUTEST(BrIndexNaiveTest, PhiPhiI)
{
    std::string s("abcdabcdabcdabcdhello");
    br_index_naive<> idx(s);
    IUTEST_EXPECT_EQ(21,idx.Phi(0));
    IUTEST_EXPECT_EQ(0,idx.Phi(4));
    IUTEST_EXPECT_EQ(4,idx.Phi(8));
    IUTEST_EXPECT_EQ(8,idx.Phi(12));
    IUTEST_EXPECT_EQ(12,idx.Phi(1));
    IUTEST_EXPECT_EQ(1,idx.Phi(5));
    IUTEST_EXPECT_EQ(5,idx.Phi(9));
    IUTEST_EXPECT_EQ(9,idx.Phi(13));
    IUTEST_EXPECT_EQ(13,idx.Phi(2));
    IUTEST_EXPECT_EQ(2,idx.Phi(6));
    IUTEST_EXPECT_EQ(6,idx.Phi(10));
    IUTEST_EXPECT_EQ(10,idx.Phi(14));
    IUTEST_EXPECT_EQ(14,idx.Phi(3));
    IUTEST_EXPECT_EQ(3,idx.Phi(7));

    IUTEST_EXPECT_EQ(0,idx.PhiI(21));
    IUTEST_EXPECT_EQ(4,idx.PhiI(0));
    IUTEST_EXPECT_EQ(8,idx.PhiI(4));
    IUTEST_EXPECT_EQ(12,idx.PhiI(8));
    IUTEST_EXPECT_EQ(1,idx.PhiI(12));
    IUTEST_EXPECT_EQ(5,idx.PhiI(1));
    IUTEST_EXPECT_EQ(9,idx.PhiI(5));
    IUTEST_EXPECT_EQ(13,idx.PhiI(9));
    IUTEST_EXPECT_EQ(2,idx.PhiI(13));
    IUTEST_EXPECT_EQ(6,idx.PhiI(2));
    IUTEST_EXPECT_EQ(10,idx.PhiI(6));
    IUTEST_EXPECT_EQ(14,idx.PhiI(10));
    IUTEST_EXPECT_EQ(3,idx.PhiI(14));
    IUTEST_EXPECT_EQ(7,idx.PhiI(3));

}



IUTEST(BrIndexNaiveTest, DNALikeTextLocate)
{
    std::string s("AAAATGCCGCCGCCATAAA");
    br_index_naive<> idx(s);

    auto vec = idx.locate("C");
    std::sort(vec.begin(),vec.end());
    print_vec<>(vec);
    IUTEST_EXPECT_EQ(6,vec[0]);
    IUTEST_EXPECT_EQ(7,vec[1]);
    IUTEST_EXPECT_EQ(9,vec[2]);
    IUTEST_EXPECT_EQ(10,vec[3]);
    IUTEST_EXPECT_EQ(12,vec[4]);
    IUTEST_EXPECT_EQ(13,vec[5]);

    vec = idx.locate("CC");
    std::sort(vec.begin(),vec.end());
    print_vec<>(vec);
    IUTEST_EXPECT_EQ(6,vec[0]);
    IUTEST_EXPECT_EQ(9,vec[1]);
    IUTEST_EXPECT_EQ(12,vec[2]);

    vec = idx.locate("CC",false);
    std::sort(vec.begin(),vec.end());
    print_vec<>(vec);
    IUTEST_EXPECT_EQ(6,vec[0]);
    IUTEST_EXPECT_EQ(9,vec[1]);
    IUTEST_EXPECT_EQ(12,vec[2]);


    vec = idx.locate("GCC");
    std::sort(vec.begin(),vec.end());
    print_vec<>(vec);
    IUTEST_EXPECT_EQ(5,vec[0]);
    IUTEST_EXPECT_EQ(8,vec[1]);
    IUTEST_EXPECT_EQ(11,vec[2]);

    vec = idx.locate("GCC",false);
    std::sort(vec.begin(),vec.end());
    print_vec<>(vec);
    IUTEST_EXPECT_EQ(5,vec[0]);
    IUTEST_EXPECT_EQ(8,vec[1]);
    IUTEST_EXPECT_EQ(11,vec[2]);
}
