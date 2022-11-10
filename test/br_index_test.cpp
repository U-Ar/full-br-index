#include "iutest.hpp"
#include <vector>
#include <fstream>
#include <string>

#include "br_index.hpp"

using namespace bri;
using namespace std;

template<class T, class U>
bool equal_set(vector<T>& v1, vector<U>& v2)
{
    if (v1.size() != v2.size()) return false;
    std::sort(v1.begin(),v1.end());
    std::sort(v2.begin(),v2.end());
    for (size_t i = 0; i < v1.size(); ++i)
    {
        if (v1[i] != v2[i]) return false;
    }
    return true;

}

IUTEST(BrIndexInmemoryTest, BasicLocate)
{
    std::string s("aaaaaaaaaaaaaaaaaaaa");
    br_index idx(s,1);
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

}

IUTEST(BrIndexInmemoryTest, PeriodicTextLocate)
{
    std::string s("abcdabcdabcdabcdhello");
    br_index idx(s,1);
    auto vec = idx.locate("abcd");
    std::sort(vec.begin(),vec.end());
    IUTEST_ASSERT_EQ(0,vec[0]);
    IUTEST_ASSERT_EQ(4,vec[1]);
    IUTEST_ASSERT_EQ(8,vec[2]);
    IUTEST_ASSERT_EQ(12,vec[3]);

}

IUTEST(BrIndexInmemoryTest, PhiPhiI)
{
    std::string s("abcdabcdabcdabcdhello");
    br_index idx(s,1);
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



IUTEST(BrIndexInmemoryTest, DNALikeTextLocate)
{
    std::string s("AAAATGCCGCCGCCATAAA");
    br_index idx(s,1);

    auto vec = idx.locate("C");
    std::sort(vec.begin(),vec.end());
    IUTEST_EXPECT_EQ(6,vec[0]);
    IUTEST_EXPECT_EQ(7,vec[1]);
    IUTEST_EXPECT_EQ(9,vec[2]);
    IUTEST_EXPECT_EQ(10,vec[3]);
    IUTEST_EXPECT_EQ(12,vec[4]);
    IUTEST_EXPECT_EQ(13,vec[5]);

    vec = idx.locate("CC");
    std::sort(vec.begin(),vec.end());
    IUTEST_EXPECT_EQ(6,vec[0]);
    IUTEST_EXPECT_EQ(9,vec[1]);
    IUTEST_EXPECT_EQ(12,vec[2]);

    vec = idx.locate("GCC");
    std::sort(vec.begin(),vec.end());
    IUTEST_EXPECT_EQ(5,vec[0]);
    IUTEST_EXPECT_EQ(8,vec[1]);
    IUTEST_EXPECT_EQ(11,vec[2]);
}

IUTEST(BrIndexInmemoryTest, ExtensionContraction)
{
    string input("abracadabra");
    br_index idx(input,1,false);
    auto init = idx.get_initial_sample();
    auto s = init;

    s = idx.left_extension(init,'a');
    {
        vector<ulint> exp{0,3,5,7,10};
        vector<ulint> loc(idx.locate_sample(s));
        IUTEST_ASSERT(equal_set<>(exp,loc));
    }
    s = idx.left_extension(s,'r');
    {
        vector<ulint> exp{2,9};
        vector<ulint> loc(idx.locate_sample(s));
        IUTEST_ASSERT(equal_set<>(exp,loc));
    }
    s = idx.left_contraction(s);
    {
        vector<ulint> exp{0,3,5,7,10};
        vector<ulint> loc(idx.locate_sample(s));
        IUTEST_ASSERT(equal_set<>(exp,loc));
    }
    s = idx.left_contraction(s);
    IUTEST_EXPECT_EQ(12,s.size());

    s = idx.right_extension(init,'r');
    {
        vector<ulint> exp{2,9};
        vector<ulint> loc(idx.locate_sample(s));
        IUTEST_ASSERT(equal_set<>(exp,loc));
    }
    s = idx.right_extension(s,'a');
    {
        vector<ulint> exp{2,9};
        vector<ulint> loc(idx.locate_sample(s));
        IUTEST_ASSERT(equal_set<>(exp,loc));
    }
    s = idx.right_contraction(s);
    {
        vector<ulint> exp{2,9};
        vector<ulint> loc(idx.locate_sample(s));
        IUTEST_ASSERT(equal_set<>(exp,loc));
    }
    s = idx.right_contraction(s);
    IUTEST_EXPECT_EQ(12,s.size());
    
}