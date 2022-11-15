#include "iutest.hpp"

#include <vector>
#include <fstream>
#include "huffman_string.hpp"

using namespace bri;

IUTEST(HuffmanStringTest,AllAText) {
    std::string s;
    for (int i = 0; i < 1000; ++i) s += 'a';
    huffman_string hf(s);

    for (int i = 0; i < 1000; ++i)
    {
        IUTEST_ASSERT_EQ(hf[i],'a');
        IUTEST_ASSERT_EQ(hf.rank(i,'a'),i);
        IUTEST_ASSERT_EQ(hf.rank(i,'b'),0);
        IUTEST_ASSERT_EQ(hf.select(i,'a'),i);
    }
    IUTEST_ASSERT_EQ(hf.rank(1000,'a'),1000);

}

IUTEST(HuffmanStringTest,PeriodicText) {
    std::string s;
    for (unsigned long i = 0; i < 255 * 100; ++i)
        s.push_back(i%255+1);
    huffman_string hf(s);

    for (uchar c = 1; c != 0; ++c)
    {
        IUTEST_ASSERT_EQ(100,hf.rank(25500,c));
        IUTEST_ASSERT_EQ(c-1,hf.select(0,c));
    }

}

IUTEST(HuffmanStringTest,Substitution) {
    std::string s;
    for (int i = 0; i < 1000; ++i) s += 'a';
    auto hf = huffman_string(s);

    for (int i = 0; i < 1000; ++i)
    {
        IUTEST_ASSERT_EQ(hf[i],'a');
        IUTEST_ASSERT_EQ(hf.rank(i,'a'),i);
        IUTEST_ASSERT_EQ(hf.rank(i,'b'),0);
        IUTEST_ASSERT_EQ(hf.select(i,'a'),i);
    }
    IUTEST_ASSERT_EQ(hf.rank(1000,'a'),1000);
}

IUTEST(HuffmanStringTest,SaveLoad) {
    std::string s;
    for (int i = 0; i < 1000; ++i) s += 'a';
    auto hf = huffman_string(s);

    std::ofstream ofs("test-tmp/huffman_string_test.tmp");
    unsigned long w_bytes = hf.serialize(ofs);
    std::cout << w_bytes << " bytes" << std::endl;
    ofs.close();

    std::ifstream ifs("test-tmp/huffman_string_test.tmp");
    hf.load(ifs);
    ifs.close();

    for (int i = 0; i < 1000; ++i)
    {
        IUTEST_ASSERT_EQ(hf[i],'a');
        IUTEST_ASSERT_EQ(hf.rank(i,'a'),i);
        IUTEST_ASSERT_EQ(hf.rank(i,'b'),0);
        IUTEST_ASSERT_EQ(hf.select(i,'a'),i);
    }
    IUTEST_ASSERT_EQ(hf.rank(1000,'a'),1000);
}