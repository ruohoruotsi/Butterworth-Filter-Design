#include "sos2tf.hpp"
#include "conv.hpp"

// implement matlab statement: b(1:3) = sos(1, 1:3);
// exclude end index.
void slice_assign(vector<double> &lhs, int lstart, int lend, vector<double> &rhs, int rstart, int rend)
{
    for (int i = 0; i < lend - lstart; i++)
    {
        lhs[lstart + i] = rhs[rstart + i];
    }
}

void sos2tf(const vector<Biquad> &sos, double gain, vector<double> &b, vector<double> &a)
{
    // Step1: Decompose b and a coefficients
    int L = sos.size();

    b.resize(2 * L + 1, 0);
    a.resize(2 * L + 1, 0);

    // initializes the first polynomial coefficient
    // b(1 : 3) = sos(1, 1 : 3);
    // a(1 : 3) = sos(1, 4 : 6);
    b[0] = sos[0].b0;
    b[1] = sos[0].b1;
    b[2] = sos[0].b2;
    a[0] = 1;
    a[1] = sos[0].a1;
    a[2] = sos[0].a2;

    // iterate every Biquad
    for (int r = 1; r < L; r++)
    {
        // get b and a of SOS
        vectord b1(3), a1(3);
        b1[0] = sos[r].b0;
        b1[1] = sos[r].b1;
        b1[2] = sos[r].b2;
        a1[0] = 1;
        a1[1] = sos[r].a1;
        a1[2] = sos[r].a2;
        // polynomial multiply by convolve.
        // b(1:2*(r+1)+1) = myconv(b(1:2*r+1), b1);
        std::vector<double>::const_iterator b_begin = b.begin();
        std::vector<double>::const_iterator b_end = b.begin() + 2 * r + 1;
        vectord b0(b_begin, b_end);
        vectord bc;
        conv(b0, b1, bc);
        slice_assign(b, 0, bc.size(), bc, 0, bc.size());

        // a(1:2*(r+1)+1) = myconv(a(1:2*r+1), a1);
        std::vector<double>::const_iterator a_begin = a.begin();
        std::vector<double>::const_iterator a_end = a.begin() + 2 * r + 1;
        vectord a0(a_begin, a_end);
        vectord ac;
        conv(a0, a1, ac);
        slice_assign(a, 0, ac.size(), ac, 0, ac.size());
    }
    // cut size to 2L+1, because the convolution adds two zeros to the end
    b.resize(2*L+1);
    a.resize(2*L+1);
    // multiply gain
    for (int i = 0; i < b.size(); i++) {
        b[i] *= gain;
    }
}