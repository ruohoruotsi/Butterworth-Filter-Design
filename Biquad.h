/*
 
 This file is part of Butterworth Filter Design, a pair C++ classes and an
 accompanying suite of unit tests for designing high order Butterworth IIR &
 EQ filters using the bilinear transform.
 The generated filter coefficients are split out into cascaded biquad sections,
 for easy use in your garden variety biquad or second-order section (SOS).
 
 Reference: http://en.wikipedia.org/wiki/Butterworth_filter
 http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 
 
 Copyright (C) 2013,  iroro orife
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 */


#include <vector>
#include <complex>

using namespace std;

// A biquad filter expression:
// y[n] = b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] + a1 * y[n-1] + a2 * y[n-2];

class Biquad {
    
public:
    Biquad();
    ~Biquad();
    
    
    // Coefficients are public, as client classes and test harnesses need direct access.
    // Special accessors to better encapsulate data, would be less readable,
    // so skipping for pedagogical reasons. Another idea would be to make this a struct.
    //
    double b0, b1, b2, a1, a2;  // second order section variable
    double a3, a4, b3, b4;      // fourth order section variables
    
    
    // Coefficients for a DF2T fourth order section (Used for EQ filters)
    void DF2TFourthOrderSection(double B0, double B1, double B2, double B3, double B4,
                                double A0, double A1, double A2, double A3, double A4);
    
    // Coefficients for a DF2T biquad section.
    void DF2TBiquad(double B0, double B1, double B2,
                    double A0, double A1, double A2);
};



class BiquadChain {
    
public:
    
    BiquadChain(int count);
    
    BiquadChain();
    ~BiquadChain();
    
    void resize(int count);
    void reset();
    
    // Process the biquad filter chain on the input buffer, write to output buffer
    // buffers arrays contain count elements with a single stride separating successive elements.
    void processBiquad(const float * input, float * output, const int stride, const int count, const Biquad * coeffs);
    
    // Extension for fourth order sections
    void processFourthOrderSections(const float * input, float * output, int stride, int count, const Biquad * coeffs);
    
private:
    
    int numFilters;
    double _xn1, _xn2;
    vector <double> _yn, _yn1, _yn2;
    
    // Fourth order section variables
    double xn3, xn4;
    vector <double> _yn3, _yn4;
    
    void allocate(int count);
};
