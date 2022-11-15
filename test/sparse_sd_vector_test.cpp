#include "iutest.hpp"

#include <vector>
#include <fstream>
#include "sparse_sd_vector.hpp"

using namespace bri;

IUTEST(SparseSdVectorTest, AllZeroBasic) {
    std::vector<bool> vec(10000,0);
    bri::sparse_sd_vector bv(vec);

    IUTEST_ASSERT_EQ(bv.size(), 10000);

    for (size_t i = 0; i < bv.size(); ++i)
    {
        IUTEST_ASSERT_EQ(0,bv[i]);
        IUTEST_ASSERT_EQ(0,bv.rank(i));
    }
    IUTEST_ASSERT_EQ(0,bv.number_of_1());
}

IUTEST(SparseSdVectorTest, AllOneBasic) {
    std::vector<bool> vec(10000,1);
    bri::sparse_sd_vector bv(vec);

    for (size_t i = 0; i < bv.size(); ++i)
    {
        IUTEST_ASSERT_EQ(1,bv[i]);
        IUTEST_ASSERT_EQ(i,bv.rank(i));
        IUTEST_ASSERT_EQ(i,bv.select(i));
    }
}

IUTEST(SparseSdVectorTest, SaveLoadBasic) {
    {
        std::vector<bool> vec(10000,false);
        bri::sparse_sd_vector bv(vec);

        std::ofstream ofs("test-tmp/sparse_sd_vector_test.tmp");
        ulint w_bytes = bv.serialize(ofs);
        std::cout << w_bytes << " bytes" << std::endl;
        ofs.close();

        std::ifstream ifs("test-tmp/sparse_sd_vector_test.tmp");
        bv.load(ifs);
        ifs.close();
        std::cout << "load complete" << std::endl;
        
        for (size_t i = 0; i < bv.size(); ++i)
        {
            IUTEST_ASSERT_EQ(bv[i],0);
            IUTEST_ASSERT_EQ(bv.rank(i),0);
        }
        IUTEST_ASSERT_EQ(bv.number_of_1(),0);
        
    }
    {
        std::vector<bool> vec(10000,true);
        bri::sparse_sd_vector bv(vec);

        std::ofstream ofs("test-tmp/sparse_sd_vector_test2.tmp");
        ulint w_bytes = bv.serialize(ofs);
        std::cout << w_bytes << " bytes" << std::endl;
        ofs.close();

        std::ifstream ifs("test-tmp/sparse_sd_vector_test2.tmp");
        bv.load(ifs);
        ifs.close();
        std::cout << "load complete" << std::endl;
        
        for (size_t i = 0; i < bv.size(); ++i)
        {
            IUTEST_ASSERT_EQ(bv[i],1);
            IUTEST_ASSERT_EQ(bv.rank(i),i);
            IUTEST_ASSERT_EQ(bv.select(i),i);
        }
        
    }
}

IUTEST(SparseSdVectorTest, SubstitutionOperator)
{
    std::vector<bool> vec(10000,1);
    auto bv = bri::sparse_sd_vector(vec);

    for (size_t i = 0; i < bv.size(); ++i)
    {
        IUTEST_ASSERT_EQ(1,bv[i]);
        IUTEST_ASSERT_EQ(i,bv.rank(i));
        IUTEST_ASSERT_EQ(i,bv.select(i));
    }
}

IUTEST(SparseSdVectorTest, Predecessor)
{
    {
        std::vector<bool> vec(10000,1);
        bri::sparse_sd_vector bv(vec);

        for (size_t i = 1 ; i < bv.size(); ++i)
        {
            IUTEST_ASSERT_EQ(bv.predecessor(i),i-1);
        }
    }
    {
        std::vector<bool> vec(10000,0);
        vec[0] = 1;
        vec[5000] = 1;
        
        bri::sparse_sd_vector bv(vec);

        for (size_t i = 1; i <= 5000; ++i)
        {
            IUTEST_ASSERT_EQ(bv.predecessor(i),0);
        }
        for (size_t i = 5001; i < 10000; ++i)
        {
            IUTEST_ASSERT_EQ(bv.predecessor(i),5000);
        }
    }
}

IUTEST(SparseSdVectorTest, PredecessorRankCircular)
{
//return rank(i) == 0 ? number_of_1() - 1 : rank(i) - 1;
    {
        std::vector<bool> vec(10000,1);
        bri::sparse_sd_vector bv(vec);

        IUTEST_ASSERT_EQ(bv.predecessor_rank_circular(0),9999);
        for (size_t i = 1 ; i < bv.size(); ++i)
        {
            IUTEST_ASSERT_EQ(bv.predecessor_rank_circular(i),i-1);
        }
    }
    {
        std::vector<bool> vec(10000,0);
        vec[0] = 1;
        vec[5000] = 1;
        
        bri::sparse_sd_vector bv(vec);

        IUTEST_ASSERT_EQ(bv.predecessor_rank_circular(0),1);
        for (size_t i = 1; i <= 5000; ++i)
        {
            IUTEST_ASSERT_EQ(bv.predecessor_rank_circular(i),0);
        }
        for (size_t i = 5001; i < 10000; ++i)
        {
            IUTEST_ASSERT_EQ(bv.predecessor_rank_circular(i),1);
        }
    }
}

IUTEST(SparseSdVectorTest, NumberOfOne)
{
    {
        std::vector<bool> vec(10000,1);
        bri::sparse_sd_vector bv(vec);
        IUTEST_ASSERT_EQ(bv.number_of_1(),10000);
    }
    {
        std::vector<bool> vec(10000,0);
        bri::sparse_sd_vector bv(vec);
        IUTEST_ASSERT_EQ(bv.number_of_1(),0);
    }
    {
        std::vector<bool> vec(10000,0);
        vec[0] = 1;
        vec[5000] = 1;
        vec[9999] = 1;
        bri::sparse_sd_vector bv(vec);
        IUTEST_ASSERT_EQ(bv.number_of_1(),3);
    }
}