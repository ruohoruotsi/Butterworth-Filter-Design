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


#include "Biquad.h"

typedef complex <double> complex_double;


class Butterworth {
    
public:
    
    Butterworth(){
    }
    ~Butterworth(){
    }
    
    static vector <complex_double> prototypeAnalogLowPass(int filterOrder);
    
    bool loPass(double fs, double f1, double f2, int filterOrder,
                vector <Biquad> & coeffs, double & overallGain){
        return coefficients(kLoPass, fs, f1, f2, filterOrder, coeffs, overallGain);
    }
    
    bool hiPass(double fs, double f1, double f2, int filterOrder,
                vector <Biquad> & coeffs, double & overallGain){
        return coefficients(kHiPass, fs, f1, f2, filterOrder, coeffs, overallGain);
    }
    
    bool bandPass(double fs, double f1, double f2, int filterOrder,
                  vector <Biquad> & coeffs, double & overallGain){
        return coefficients(kBandPass, fs, f1, f2, filterOrder, coeffs, overallGain);
    }
    
    bool bandStop(double fs, double f1, double f2, int filterOrder,
                  vector <Biquad> & coeffs, double & overallGain){
        return coefficients(kBandStop, fs, f1, f2, filterOrder, coeffs, overallGain);
    }
    
    
    enum FILTER_TYPE {
        kLoPass     = 10000,
        kHiPass     = 10001,
        kBandPass   = 10002,
        kBandStop   = 10003,
        kLoShelf    = 10004,
        kHiShelf    = 10005, // high order EQ
        kParametric = 10006  // high order EQ
    };
    
    
    //******************************************************************************
    
    // generic coeffs
    bool coefficients(FILTER_TYPE filter, const double fs, const double freq1_cutoff, const double freq2_cutoff, const int filterOrder,
                      vector <Biquad> & coeffs, double & overallGain);
    
    // high order EQ
    bool coefficientsEQ(FILTER_TYPE filter, const double fs, const double f1, const double f2, const int filterOrder,
                        vector <Biquad> & coeffs, double overallGain);
    
    
private:
    
    double blt(complex_double & sz);
    bool s2Z();
    bool zp2SOS();
    
    // Private helper methods to convert the lowpass analog prototype to the desired filter type
    void convert2lopass();
    void convert2hipass();
    void convert2bandpass();
    void convert2bandstop();
    
    // Internal state used during computation of coefficients
    double f1, f2;
    
    int numPoles, numZeros;
    vector <complex_double> zeros;
    vector <complex_double> poles;
    
    double Wc; // Omega cutoff == passband edge freq
    double bw; // Bandwidth
    
    double gain;
    double preBLTgain;
    
    int nba;
    double * ba;
};

