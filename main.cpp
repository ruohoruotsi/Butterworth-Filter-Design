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


#include <iostream>
#include <sndfile.hh>
#include "Butterworth.h"

#define CATCH_CONFIG_MAIN   // Tells Catch to provide a main()
#include "catch.hpp"        // Catch Unit test framework

using namespace std;


TEST_CASE("Analogue lowpass prototype generation", " Verify poles match MATLAB buttap(8)"){
    
    //******************************************************************************
    //
    // baseline_poles values generated with MATLAB command >> buttap(8);
    //
    //******************************************************************************
    
    int filterOrder = 8;
    
    vector <complex_double> poles;
    poles.resize(filterOrder);
    vector <complex_double> baseline_poles;
    baseline_poles.resize(filterOrder);
    
    // Butterworth butter;
    vector <complex_double> tempPoleStorage = Butterworth::prototypeAnalogLowPass(filterOrder);
    
    // copy tmppole into poles
    int index = 0;
    for(vector <complex_double>::iterator itr = tempPoleStorage.begin(); itr != tempPoleStorage.end(); itr++){
        poles[index] = *itr;
        index++;
    }
    
    baseline_poles[0] = complex_double(-0.19509032201613,  0.98078528040323);
    baseline_poles[1] = complex_double(-0.19509032201613, -0.98078528040323);
    baseline_poles[2] = complex_double(-0.55557023301960,  0.83146961230255);
    baseline_poles[3] = complex_double(-0.55557023301960, -0.83146961230255);
    baseline_poles[4] = complex_double(-0.83146961230255,  0.55557023301960);
    baseline_poles[5] = complex_double(-0.83146961230255, -0.55557023301960);
    baseline_poles[6] = complex_double(-0.98078528040323,  0.19509032201613);
    baseline_poles[7] = complex_double(-0.98078528040323, -0.19509032201613);
    
    // Tolerance
    const double EPSILON = 1.0e-14;
    
    // Test
    for(uint32_t i = 0; i < filterOrder; i++){
        CAPTURE(i);
        REQUIRE(abs(baseline_poles[i].real() - poles[i].real()) <= EPSILON);
    }
}


TEST_CASE("Lowpass filter coefficients", " Verify coefficients for a lowpass filter"){
    
    //******************************************************************************
    // MATLAB (vR14) coefficients generated with the following code:
    //
    // [z, p, k] = butter(8, 500, 's');
    // [Zd, Pd, Kd] = bilinear(z, p, k, 44100);
    // [sos, g] = zpk2sos(Zd, Pd, Kd)
    //******************************************************************************
    
    int filterOrder = 8;
    double overallGain = 1.0;
    
    vector <Biquad> coeffs;  // second-order sections (sos)
    Butterworth butterworth;
    
    bool designedCorrectly = butterworth.loPass(44100,  // fs
                                                500,    // freq1
                                                0,      // freq2. N/A for lowpass
                                                filterOrder,
                                                coeffs,
                                                overallGain);
    REQUIRE(designedCorrectly == true);
    
    //******************************************************************************
    // MATLAB coefficients: first section
    //******************************************************************************
    
    double b0 = 1.0;
    double b1 = 2.0;
    double b2 = 1.0;
    // double a0 = 1.0;
    double a1 = -1.96762058043629 * (-1); // to convert to DF2T (direct form II transpose)
    double a2 = 0.97261960500367  * (-1); // to convert to DF2T (direct form II transpose)
    
    int i = 0;
    
    const double EPSILON = 1.0e-4;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    //******************************************************************************
    //  MATLAB coefficients: second section
    //******************************************************************************
    
    b0 = 1.0;
    b1 = 2.0;
    b2 = 1.0;
    // a0 = 1.0;
    a1 = -1.91907529383978 * (-1);
    a2 = 0.92395098208778  * (-1);
    
    i = 1;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    //******************************************************************************
    //  MATLAB coefficients: third section
    //******************************************************************************
    
    b0 = 1.0;
    b1 = 2.0;
    b2 = 1.0;
    // a0 = 1.0;
    a1 = -1.88350864178159 * (-1);
    a2 = 0.88829396780773  * (-1);
    
    i = 2;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    //******************************************************************************
    //  MATLAB coefficients: fourth section
    //******************************************************************************
    
    b0 = 1.0;
    b1 = 2.0;
    b2 = 1.0;
    // a0 = 1.0;
    a1 = -1.86480445083537 * (-1);
    a2 = 0.86954225616013  * (-1);
    
    i = 3;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
}


TEST_CASE("Highpass filter coefficients", " Verify coefficients for a highpass filter"){
    
    //******************************************************************************
    // MATLAB (vR14) coefficients generated with the following code:
    //
    // [z, p, k] = butter(8, 500, 'high', 's');
    // [Zd, Pd, Kd] = bilinear(z, p, k, 44100);
    // [sos, g] = zpk2sos(Zd, Pd, Kd)
    //******************************************************************************
    
    int filterOrder = 8;
    double overallGain = 1.0;
    const double EPSILON = 1.0e-4;
    
    vector <Biquad> coeffs;  // second-order sections (sos)
    Butterworth butterworth;
    
    bool designedCorrectly = butterworth.hiPass(44100,  // fs
                                                500,    // freq1
                                                0,      // freq2
                                                filterOrder,
                                                coeffs,
                                                overallGain);
    REQUIRE(designedCorrectly == true);
    
    //******************************************************************************
    //  MATLAB coefficients: first section
    //******************************************************************************
    
    double b0 = 1.0;
    double b1 = -2.0;
    double b2 = 1.0;
    // double a0 = 1.0;
    double a1 = -1.96762058043629 * (-1);
    double a2 = 0.97261960500367  * (-1);
    
    int i = 0;
    
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    //******************************************************************************
    //  MATLAB coefficients: second section
    //******************************************************************************
    
    b0 = 1.0;
    b1 = -2.0;
    b2 = 1.0;
    // a0 = 1.0;
    a1 = -1.91907529383978 * (-1);
    a2 = 0.92395098208778  * (-1);
    
    i = 1;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    //******************************************************************************
    //  MATLAB coefficients: third section
    //******************************************************************************
    
    b0 = 1.0;
    b1 = -2.0;
    b2 = 1.0;
    // a0 = 1.0;
    a1 = -1.88350864178159 * (-1);
    a2 = 0.88829396780773  * (-1);
    
    i = 2;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    //******************************************************************************
    //  MATLAB coefficients:  fourth section
    //******************************************************************************
    
    b0 = 1.0;
    b1 = -2.0;
    b2 = 1.0;
    // a0 = 1.0;
    a1 = -1.86480445083537 * (-1);
    a2 = 0.86954225616013  * (-1);
    
    i = 3;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
}


TEST_CASE("Parametric EQ coefficients", " Verify coefficients for a parametric filter"){
    
    //******************************************************************************
    // MATLAB coefficients were generated with the following code in
    // MATLAB R14 with code hpeq1.m from
    //
    // Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer
    //            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
    //******************************************************************************
    
    int filterOrder = 8;
    double overallGain = 12.0;
    const double EPSILON = 1.0e-4;
    
    vector <Biquad> coeffs;
    Butterworth butterworth;
    
    bool designedCorrectly = butterworth.coefficientsEQ(Butterworth::kParametric,
                                                        40000,      // fs
                                                        3500,       // freq1
                                                        6500,       // freq2   Df = 3kHz
                                                        filterOrder,
                                                        coeffs,
                                                        overallGain);
    
    REQUIRE(designedCorrectly == true);
    
    //******************************************************************************
    //  MATLAB coefficients: first section
    //******************************************************************************
    
    double b0 = 1.0356;
    double b1 = -2.6699;
    double b2 = 3.4403;
    double b3 = -2.3905;
    double b4 = 0.8435;
    // double a0 = 1.0;
    double a1 = -2.6478 * (-1);
    double a2 = 3.4810  * (-1);
    double a3 = -2.4127 * (-1);
    double a4 = 0.8384  * (-1);
    
    int i = 0;
    
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(b3 - coeffs[i].b3) <= EPSILON);
    REQUIRE(abs(b4 - coeffs[i].b4) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    REQUIRE(abs(a3 - coeffs[i].a3) <= EPSILON);
    REQUIRE(abs(a4 - coeffs[i].a4) <= EPSILON);
    
    //******************************************************************************
    //  MATLAB coefficients: second section
    //******************************************************************************
    
    b0 = 1.0555;
    b1 = -2.5476;
    b2 = 2.9933;
    b3 = -1.8553;
    b4 = 0.5794;
    //  a0 = 1.0;
    a1 = -2.4927 * (-1);
    a2 = 3.0287  * (-1);
    a3 = -1.9102 * (-1);
    a4 = 0.5995  * (-1);
    
    i = 1;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(b3 - coeffs[i].b3) <= EPSILON);
    REQUIRE(abs(b4 - coeffs[i].b4) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    REQUIRE(abs(a3 - coeffs[i].a3) <= EPSILON);
    REQUIRE(abs(a4 - coeffs[i].a4) <= EPSILON);
}


TEST_CASE("Low shelf EQ coefficients", " Verify coefficients for a low shelf filter"){
    
    //******************************************************************************
    // MATLAB coefficients were generated with the following code in
    // MATLAB R14 with code hpeq1.m from
    //
    // Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer
    //            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
    //******************************************************************************
    
    int filterOrder = 4;
    double overallGain = 12.0;
    const double EPSILON = 1.0e-4;
    
    vector <Biquad> coeffs;  // second-order sections (sos)
    Butterworth butterworth;
    bool designedCorrectly = butterworth.coefficientsEQ(Butterworth::kLoShelf,
                                                        40000,      // fs
                                                        0,          // freq1
                                                        13000,      // freq2   Df = 13kHz
                                                        filterOrder,
                                                        coeffs,
                                                        overallGain);
    
    REQUIRE(designedCorrectly == true);
    
    //******************************************************************************
    //  MATLAB coefficients: first section
    //******************************************************************************
    
    double b0 = 1.6392;
    double b1 = 1.7241;
    double b2 = 0.9170;
    // double a0 = 1.0;
    double a1 = 0.6565 * (-1);
    double a2 = 0.4887 * (-1);
    
    int i = 0;
    
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    //******************************************************************************
    //  MATLAB coefficients: second section
    //******************************************************************************
    
    b0 = 1.5790;
    b1 = 1.2663;
    b2 = 0.2984;
    // a0 = 1.0;
    a1 = 0.4822 * (-1);
    a2 = 0.0934 * (-1);
    
    i = 1;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    
    //******************************************************************************
    // test another set of lowshelf design requirements
    //******************************************************************************
    
    designedCorrectly = butterworth.coefficientsEQ(Butterworth::kLoShelf,
                                                   40000,   // fs
                                                   2000,    // freq1
                                                   2800,    // freq2   Df = 0.8kHz
                                                   filterOrder,
                                                   coeffs,
                                                   overallGain);
    
    REQUIRE(designedCorrectly == true);
    
    //******************************************************************************
    //  MATLAB coefficients: first section
    //******************************************************************************
    
    b0 = 1.0222;
    b1 = -1.8880;
    b2 = 0.8949;
    // a0 = 1.0;
    a1 = -1.8953 * (-1);
    a2 = 0.9099  * (-1);
    
    i = 0;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    
    //******************************************************************************
    //  MATLAB coefficients: second section
    //******************************************************************************
    
    b0 = 1.0456;
    b1 = -1.7749;
    b2 = 0.7566;
    // a0 = 1.0;
    a1 = -1.7817 * (-1);
    a2 = 0.7954 * (-1);
    
    i = 1;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
}


//-------------------------------------------------------------------------------------------------
TEST_CASE("High shelf EQ coefficients", " Verify coefficients for a high shelf filter"){
    
    //******************************************************************************
    // MATLAB coefficients were generated with the following code in
    // MATLAB R14 with code hpeq1.m from
    //
    // Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer
    //            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
    //******************************************************************************
    
    int filterOrder = 4;
    double overallGain = 12.0;
    const double EPSILON = 1.0e-4;
    
    
    vector <Biquad> coeffs;  // second-order sections (sos)
    Butterworth butterworth;
    bool designedCorrectly = butterworth.coefficientsEQ(Butterworth::kHiShelf,
                                                        40000,      // fs
                                                        0,          // freq1
                                                        13000,      // freq2   Df = 13kHz
                                                        filterOrder,
                                                        coeffs,
                                                        overallGain);
    
    
    REQUIRE(designedCorrectly == true);
    
    //******************************************************************************
    //  MATLAB coefficients: first section
    //******************************************************************************
    
    double b0 = 1.6392;
    double b1 = -1.7241;
    double b2 = 0.9170;
    // double a0 = 1.0;
    double a1 = -0.6565 * (-1);
    double a2 = 0.4887 * (-1);
    
    int i = 0;
    
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    //******************************************************************************
    //  MATLAB coefficients: second section
    //******************************************************************************
    
    b0 = 1.5790;
    b1 = -1.2663;
    b2 = 0.2984;
    // a0 = 1.0;
    a1 = -0.4822 * (-1);
    a2 = 0.0934  * (-1);
    
    i = 1;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    
    
    //******************************************************************************
    // test another set of highshelf design requirements
    //******************************************************************************
    
    designedCorrectly = butterworth.coefficientsEQ(Butterworth::kHiShelf,
                                                   40000,   // fs
                                                   2000,    // freq1
                                                   2800,    // freq2   Df = 13kHz
                                                   filterOrder,
                                                   coeffs,
                                                   overallGain);
    
    REQUIRE(designedCorrectly == true);
    
    //******************************************************************************
    //  MATLAB coefficients: first section
    //******************************************************************************
    
    b0 = 1.0222;
    b1 = 1.8880;
    b2 = 0.8949;
    // a0 = 1.0;
    a1 = 1.8953 * (-1);
    a2 = 0.9099 * (-1);
    
    i = 0;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
    
    //******************************************************************************
    //  MATLAB coefficients: second section
    //******************************************************************************
    
    b0 = 1.0456;
    b1 = 1.7749;
    b2 = 0.7566;
    // a0 = 1.0;
    a1 = 1.7817 * (-1);
    a2 = 0.7954 * (-1);
    
    i = 1;
    REQUIRE(abs(b0 - coeffs[i].b0) <= EPSILON);
    REQUIRE(abs(b1 - coeffs[i].b1) <= EPSILON);
    REQUIRE(abs(b2 - coeffs[i].b2) <= EPSILON);
    REQUIRE(abs(a1 - coeffs[i].a1) <= EPSILON);
    REQUIRE(abs(a2 - coeffs[i].a2) <= EPSILON);
}


//-------------------------------------------------------------------------------------------------
#define	BUFFER_LEN 1024
#define	MAX_CHANNELS 6

static void process_data (double *data, int count, int channels)
{
    int k;
    
    /* Process the data here.
     ** If the soundfile contains more then 1 channel you need to take care of
     ** the data interleaving youself.
     ** Current we just apply a channel dependant gain.
     */
    
    for (k = 0 ; k < count ; k+= 1)
        data [k] *= 1.0;
    
    return ;
}

static void read_write_file (const char	*infilename, const char	*outfilename)
{
    static double data [BUFFER_LEN] ;
    SNDFILE	*infile, *outfile ;
    SF_INFO		sfinfo ;
    int			readcount ;
    memset (&sfinfo, 0, sizeof (sfinfo)) ;
    
    if (! (infile = sf_open (infilename, SFM_READ, &sfinfo)))
    {
        printf ("Not able to open input file %s.\n", infilename) ;
        puts (sf_strerror (NULL)) ;
        return;
    }
    
    if (sfinfo.channels > MAX_CHANNELS)
    {	printf ("Not able to process more than %d channels\n", MAX_CHANNELS) ;
        return;
    }
    if (! (outfile = sf_open (outfilename, SFM_WRITE, &sfinfo)))
    {	printf ("Not able to open output file %s.\n", outfilename) ;
        puts (sf_strerror (NULL)) ;
        return;
    }
    
    while ((readcount = sf_read_double (infile, data, BUFFER_LEN)))
    {	process_data (data, readcount, sfinfo.channels) ;
        sf_write_double (outfile, data, readcount) ;
    }
    
    sf_close (infile); sf_close (outfile) ;
    return;
}


TEST_CASE("Test gain control issue #3"){
    
    //******************************************************************************
    // For this test you will need libsndfile.a which lives in
    // /opt/local/lib/libsndfile.a installed on MacOS via macports
    // https://trac.macports.org/browser/trunk/dports/audio/libsndfile/Portfile
    //******************************************************************************
    
    const char * infilename = "sweep_0_20K.wav" ;
    const char * outfilename = "sweep_0_20K_filtered.wav" ;
    read_write_file(infilename,outfilename);
}
