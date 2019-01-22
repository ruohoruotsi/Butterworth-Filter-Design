## Butterworth Filter Design

Butterworth Filter Design is a collection of C++ classes and an accompanying suite of unit tests for designing high order Butterworth IIR and EQ filters using the bilinear transform. 

The generated filter coefficients are split out into cascaded biquad sections, for easy use in any biquad or second-order section (SOS) implementation.

### Features
 * Lowpass, Highpass, Bandpass and Bandstop IIR & EQ filter design 
 * Low and High Shelving filter design
 * High order Parametric boost/cut EQ filter design
 *Biquad and Biquad Chain implementations (for filtering audio buffers with cascaded biquad sections)
 * Compact, readable and well-commented codebase


### Unit tests
As with any good audio signal processing toolkit, there are unit tests that provide basic proof of correctness. There are currently 6 primary test cases that check 113 different conditions.

Unit tests live in `main.cpp` and are written using the compact [Catch](https://github.com/philsquared/Catch) test framework for C++.

### Prerequisites

 * [SCONS](http://scons.org) as a cross-platform build system to build, test and run examples. 
    * On MacOS use [Homebrew](https://brew.sh): `$brew install scons` or [MacPorts](https://www.macports.org) `port install scons`
    * On Linux: `apt-get install scons`
 * [libsndfile](http://www.mega-nerd.com/libsndfile): `brew install libsndfile`

### Usage
The unit tests are a good place to start for a survey of usage. 

For example, to design an 8-tap Butterworth lowpass filter with a cutoff @ 500Hz (which will generate coefficients for 4 biquad filters), running @ 44100Hz, with unity gain (1.0) execute the following:


```
vector <Biquad> coeffs;  // array of biquad filters (for this case, array size = 4 )
Butterworth butterworth;
bool designedCorrectly = butterworth.loPass(44100,  // fs
					    500,    // freq1
					    0,      // freq2. N/A for lowpass
					    8, 	    // filter order,
					    coeffs, // coefficient array being filled
					    1.0);   // overall gain
```



To generate the same set of coefficients in MATLAB (R14) as a comparison, to double-check our work, execute the following MATLAB commands:

```
[z, p, k] = butter(8, 500, 's');		% designs a 8-tap lowpass s-domain filter
[Zd, Pd, Kd] = bilinear(z, p, k, 44100);	% analog-to-digital filter conversion
[sos, g] = zpk2sos(Zd, Pd, Kd)			% zero-pole-gain form to second-order sections (SOS)
```

												

### Other filter design repos on GitHub
* Vinnie Falco	[DSP Filters](https://github.com/vinniefalco/DSPFilters)
* Mofr	[Butterworth](https://github.com/mofr/Butterworth)

### Terms and Conditions
For my [Master's thesis](https://github.com/ruohoruotsi/Riddim), and in the course of work writing audio plugins and music apps, I've implemented a few different IIR and EQ filter design classes (e.g. [RBJ's EQ Cookbook](http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt), Cheby and Butter IIRs). These classes are a fresh rewrite and should be considered unoptimized reference code, with an emphasis on clarity (and for pedagogical reasons, being able to refer to the textbook, see the math and understand the code). For this reason, in the hope they can be useful, these classes are provided under [GPL v3](http://www.gnu.org/licenses/gpl.html). 

Butterworth Filter Design - Copyright © 2013 &nbsp; [iroro orife](http://github.com/ruohoruotsi)

Source code is provided under [GPL v3](http://www.gnu.org/licenses/gpl.html)
