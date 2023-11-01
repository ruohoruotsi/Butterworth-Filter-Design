#ifndef SOS2TF_H
#define SOS2TF_H

#include "Biquad.h"

/**
 * My 2nd-order sections to transfer function model conversion.
 * See matlab `sos2tf`.
 * 
 *   [B,A] = MYSOS2TF(SOS,GAIN) returns the numerator and denominator
 *   coefficients B and A of the discrete-time linear system given
 *   by the gain G and the matrix SOS in second-order sections form.
 *
 *   SOS is an L by 6 matrix which contains the coefficients of each
 *   second-order section in its rows:
 *       SOS = [ b01 b11 b21  1 a11 a21
 *               b02 b12 b22  1 a12 a22
 *               ...
 *               b0L b1L b2L  1 a1L a2L ]
 *
 *   The system transfer function is the product of the second-order
 *   transfer functions of the sections times the gain G. If G is not
 *   specified, it defaults to 1. Each row of the SOS matrix describes
 *   a 2nd order transfer function as
 *       b0k +  b1k z^-1 +  b2k  z^-2
 *       ----------------------------
 *       1 +  a1k z^-1 +  a2k  z^-2
 *   where k is the row index.
 *
 * In each row of the sos matrix, the first three elements form the
 * coefficients of the numerator polynomial of a transfer function and
 * the next three elements form the coefficients of the denominator
 * polynomial. B and A are formed by multiplying all numerator and
 * denominator polynomials respectively. Since each polynomial is of
 * degree two, multiplying L such polynomials results in a polynomial of
 * degree 2*L.
 *
 * Use convolution to multiply two polynomials.
 * w = conv(u,v) Returns the convolution of the vectors u and v. 
 * If u and v are vectors of polynomial coefficients, convolving them is equivalent to multiplying the two polynomials.
 *
 */
void sos2tf(const vector<Biquad> &sos, double gain, vector<double> &b, vector<double> &a);

#endif