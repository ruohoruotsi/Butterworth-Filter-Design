#ifndef CONV_H
#define CONV_H

#include <vector>

typedef std::vector<double> vectord;

/**
 * Convolution of u and v vector.
 * Return result vector c.
 * length(c) = length(u) + length(v) - 1
 */
void conv(const vectord& u, const vectord& v, vectord& c);

#endif