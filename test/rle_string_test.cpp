#include "iutest.hpp"
#include <vector>
#include <fstream>
#include <string>

#include "rle_string.hpp"

using namespace bri;

IUTEST(RleStringTest, BasicStats)
{
    std::string s;
    for (int i = 0; i < 1000; ++i) s.push_back('a');
    rle_string<> rl(s);

    IUTEST_ASSERT_EQ(1000,rl.size());
    IUTEST_ASSERT_EQ(1,rl.number_of_runs());
    ulint bytes = rl.print_space();
    ulint bytes2 = rl.get_space();
    std::cout << "RleString total: " << bytes << " bytes" << std::endl;
    IUTEST_ASSERT_EQ(bytes, bytes2);
}

IUTEST(RleStringTest, BasicRunAt)
{
    std::string s;
    for (int i = 0; i < 1000; ++i) s.push_back('a');
    rle_string<> rl(s);

    IUTEST_ASSERT_EQ(1000,rl.run_at(0));
}

IUTEST(RleStringTest, BasicRunOf)
{
    std::string s;
    for (int i = 0; i < 1000; ++i) s.push_back('a');
    rle_string<> rl(s);

    IUTEST_ASSERT_EQ(1000,rl.size());
    for (ulint i = 0; i < rl.size(); ++i)
    {
        range_t rg = rl.run_of(i);
        IUTEST_ASSERT_EQ(0,rg.first);
        IUTEST_ASSERT_EQ(rl.size()-1,rg.second);
    }
}

IUTEST(RleStringTest, BasicIndexOperator)
{
    std::string s;
    for (int i = 0; i < 1000; ++i) s.push_back('a');
    rle_string<> rl(s);

    for (ulint i = 0; i < rl.size(); ++i)
    {
        IUTEST_ASSERT_EQ('a',rl[i]);
    }
}

IUTEST(RleStringTest, BasicRank)
{
    std::string s;
    for (int i = 0; i < 1000; ++i) s.push_back('a');
    rle_string<> rl(s);

    for (ulint i = 0; i < rl.size(); ++i)
    {
        IUTEST_ASSERT_EQ(i, rl.rank(i,'a'));
        IUTEST_ASSERT_EQ(0, rl.rank(i,'b'));
    }
}

IUTEST(RleStringTest, BasicSelect)
{
    std::string s;
    for (int i = 0; i < 1000; ++i) s.push_back('a');
    rle_string<> rl(s);

    IUTEST_ASSERT_EQ(1000,rl.size());
    for (ulint i = 0; i < rl.size(); ++i)
    {
        IUTEST_ASSERT_EQ(i, rl.select(i,'a'));
    }
}

IUTEST(RleStringTest, BasicRunOfPosition)
{
    std::string s;
    for (int i = 0; i < 1000; ++i) s.push_back('a');
    rle_string<> rl(s);

    for (ulint i = 0; i < rl.size(); ++i)
    {
        IUTEST_ASSERT_EQ(0, rl.run_of_position(i));
    }
    range_t r1 = rl.run_range(0);
    IUTEST_ASSERT_EQ(0,r1.first);
    IUTEST_ASSERT_EQ(rl.size()-1,r1.second);
}

IUTEST(RleStringTest, BasicSaveLoad)
{

    std::string s;
    for (int i = 0; i < 1000; ++i) s.push_back('a');
    rle_string<> rltmp(s);

    std::ofstream ofs("test-tmp/rle_string_test.tmp");
    rltmp.serialize(ofs);
    ofs.close();

    rle_string<> rl;
    std::ifstream ifs("test-tmp/rle_string_test.tmp");
    rl.load(ifs);
    ifs.close();

    IUTEST_ASSERT_EQ(1000,rl.size());
    IUTEST_ASSERT_EQ(1,rl.number_of_runs());
    ulint bytes = rl.print_space();
    ulint bytes2 = rl.get_space();
    std::cout << "RleString total: " << bytes << " bytes" << std::endl;
    IUTEST_ASSERT_EQ(bytes, bytes2);

    IUTEST_ASSERT_EQ(1000,rl.run_at(0));

    for (ulint i = 0; i < rl.size(); ++i)
    {
        range_t rg = rl.run_of(i);
        IUTEST_ASSERT_EQ(0,rg.first);
        IUTEST_ASSERT_EQ(rl.size()-1,rg.second);
        IUTEST_ASSERT_EQ('a',rl[i]);
        IUTEST_ASSERT_EQ(i, rl.rank(i,'a'));
        IUTEST_ASSERT_EQ(0, rl.rank(i,'b'));
        IUTEST_ASSERT_EQ(i, rl.select(i,'a'));
        IUTEST_ASSERT_EQ(0, rl.run_of_position(i));
    }
    range_t r1 = rl.run_range(0);
    IUTEST_ASSERT_EQ(0,r1.first);
    IUTEST_ASSERT_EQ(rl.size()-1,r1.second);

}

IUTEST(RleStringTest, BasicSubstitution)
{

    std::string s;
    for (int i = 0; i < 1000; ++i) s.push_back('a');
    rle_string<> rltmp(s);

    rle_string<> rl;
    rl = rltmp;

    IUTEST_ASSERT_EQ(1000,rl.size());
    IUTEST_ASSERT_EQ(1,rl.number_of_runs());
    ulint bytes = rl.print_space();
    ulint bytes2 = rl.get_space();
    std::cout << "RleString total: " << bytes << " bytes" << std::endl;
    IUTEST_ASSERT_EQ(bytes, bytes2);

    IUTEST_ASSERT_EQ(1000,rl.run_at(0));

    for (ulint i = 0; i < rl.size(); ++i)
    {
        range_t rg = rl.run_of(i);
        IUTEST_ASSERT_EQ(0,rg.first);
        IUTEST_ASSERT_EQ(rl.size()-1,rg.second);
        IUTEST_ASSERT_EQ('a',rl[i]);
        IUTEST_ASSERT_EQ(i, rl.rank(i,'a'));
        IUTEST_ASSERT_EQ(0, rl.rank(i,'b'));
        IUTEST_ASSERT_EQ(i, rl.select(i,'a'));
        IUTEST_ASSERT_EQ(0, rl.run_of_position(i));
    }
    range_t r1 = rl.run_range(0);
    IUTEST_ASSERT_EQ(0,r1.first);
    IUTEST_ASSERT_EQ(rl.size()-1,r1.second);

}

IUTEST(RleStringTest, AllCharTextOperations)
{
    std::string s;
    for (ulint i = 1; i < 256; ++i) s.push_back((uchar)i);
    rle_string<> rl(s);

    IUTEST_ASSERT_EQ(255,rl.size());
    IUTEST_ASSERT_EQ(255,rl.number_of_runs());
    ulint bytes = rl.print_space();
    ulint bytes2 = rl.get_space();
    std::cout << "RleString total: " << bytes << " bytes" << std::endl;
    IUTEST_ASSERT_EQ(bytes, bytes2);

    for (ulint i = 0; i < rl.size(); ++i)
    {
        IUTEST_ASSERT_EQ(1,rl.run_at(i));

        range_t r1 = rl.run_range(i);
        IUTEST_ASSERT_EQ(i,r1.first);
        IUTEST_ASSERT_EQ(i,r1.second);

        range_t rg = rl.run_of(i);
        IUTEST_ASSERT_EQ(i,rg.first);
        IUTEST_ASSERT_EQ(i,rg.second);

        IUTEST_ASSERT_EQ(uchar(i+1),rl[i]);
        IUTEST_ASSERT_EQ(i, rl.run_of_position(i));
    }
    for (ulint i = 1; i < 256; ++i)
    {
        IUTEST_ASSERT_EQ(1, rl.rank(rl.size(),(uchar)i));
        IUTEST_ASSERT_EQ(i-1, rl.select(0,(uchar)i));
    }
}

IUTEST(RleStringTest, TwoCharTextOperations)
{
    std::string s;
    // length 10000
    for (ulint i = 0; i < 50; ++i)
    {
        for (ulint j = 0; j < 100; ++j)
            s.push_back('a');
        for (ulint j = 0; j < 100; ++j)
            s.push_back('b');
    }
    rle_string<> rl(s);

    IUTEST_ASSERT_EQ(10000,rl.size());
    IUTEST_ASSERT_EQ(100,rl.number_of_runs());
    ulint bytes = rl.print_space();
    ulint bytes2 = rl.get_space();
    std::cout << "RleString total: " << bytes << " bytes" << std::endl;
    IUTEST_ASSERT_EQ(bytes, bytes2);

    for (ulint i = 0; i < rl.size(); ++i)
    {

        range_t rg = rl.run_of(i);
        IUTEST_ASSERT_EQ(i/100,rg.first);
        IUTEST_ASSERT_EQ(i/100 * 100 + 99,rg.second);

        if (i % 200 < 100)
        {
            IUTEST_ASSERT_EQ(uchar('a'),rl[i]);

            IUTEST_ASSERT_EQ(i/200*100+i%200,rl.rank(i,'a'));
            IUTEST_ASSERT_EQ(i,rl.select(rl.rank(i,'a'),'a'));
        }
        else 
        {  
            IUTEST_ASSERT_EQ(uchar('b'),rl[i]);

            IUTEST_ASSERT_EQ(i/200*100+i%200-100,rl.rank(i,'b'));
            IUTEST_ASSERT_EQ(i,rl.select(rl.rank(i,'b'),'b'));
        }
        IUTEST_ASSERT_EQ(i/100, rl.run_of_position(i));
    }

    for (ulint i = 0; i < 100; ++i)
    {
        
        IUTEST_ASSERT_EQ(100,rl.run_at(i));

        range_t r1 = rl.run_range(i);
        IUTEST_ASSERT_EQ(i * 100,r1.first);
        IUTEST_ASSERT_EQ(i * 100 + 99,r1.second);

    }

}