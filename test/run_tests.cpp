#include "iutest.hpp"

int main(int ac, char** av)
{
    IUTEST_INIT(&ac,av);
    return IUTEST_RUN_ALL_TESTS();
}