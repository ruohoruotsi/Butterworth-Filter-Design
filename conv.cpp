
#include "conv.hpp"

/**
 * Convolution of u and v vector.
 */
void conv(const vectord &u, const vectord &v, vectord &c)
{
    // fill with zeros
    c.resize(u.size() + v.size() - 1, 0);
    // use convolution machine implements convolution
    for (int n = 0; n < c.size(); n++)
    {
        // iterate input signal u, kernal is v
        for (int k = 0; k < v.size(); k++)
        {
            // just as the Math Equation of Convolution
            int iu = n - k;
            double uu;
            if (iu < 0 || iu >= u.size())
            {
                uu = 0;
            }
            else
            {
                uu = u[iu];
            }
            // iterate kernel and cumulate it
            c[n] = c[n] + v[k] * uu;
        }
    }
}
