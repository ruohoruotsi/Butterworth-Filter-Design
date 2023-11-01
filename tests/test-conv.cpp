#define CATCH_CONFIG_MAIN // Tells Catch to provide a main()
#include "catch.hpp"      // Catch Unit test framework
#include <cstring>
#include "conv.hpp"

using namespace std;

const double EPSILON = 1E-4;

TEST_CASE("Convolution u and v", " Verify convolution is correct")
{
    vectord u = {1, 2, 3, 4, 5};
    vectord v = {1, 2, 3};
    vectord c;
    conv(u, v, c);

    vectord r = {1, 4, 10, 16, 22, 22, 15};
    for (size_t i = 0; i < c.size(); i++)
    {
        REQUIRE(abs(r[i] - c[i]) <= EPSILON);
    }
}