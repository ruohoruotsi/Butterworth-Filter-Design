#define CATCH_CONFIG_MAIN // Tells Catch to provide a main()
#include "catch.hpp"      // Catch Unit test framework
#include <cstring>
#include "conv.hpp"
#include "sos2tf.hpp"

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

TEST_CASE("SOS To Transfer Function coefficients", " Verify coefficients is correct")
{
    vector<Biquad> sos(4);
    sos[0] = Biquad(1, 2, 1, -1.11879834477625, 0.336688180692519);
    sos[1] = Biquad(1, 2, 1, -1.35028339447298, 0.660543824885074);
    sos[2] = Biquad(1, -2, 1, -1.93859735792819, 0.939736057085694);
    sos[3] = Biquad(1, -2, 1, -1.97691413926085, 0.977914024909786);

    double gain = 0.00404923315402386;

    vectord b, a;
    sos2tf(sos, gain, b, a);

    // matlab b,a
    vectord mb = {0.00404923315402386, 0, -0.0161969326160954, 0, 0.0242953989241431, 0, -0.0161969326160954, 0, 0.00404923315402386};
    vectord ma = {1, -6.38459323643827, 17.9257353790401, -28.9644589090336, 29.5037393371934, -19.4170269899449, 8.06393935848502, -1.93171376981280, 0.204378907481348};

    for (int i = 0; i < b.size(); i++)
    {
        CAPTURE(i)
        REQUIRE(abs(b[i] - mb[i]) <= EPSILON);
        REQUIRE(abs(a[i] - ma[i]) <= EPSILON);
    }
}
