## Butterworth Filter Design
=========================

Butterworth Filter Design is a collection of C++ classes and an accompanying suite of unit tests for designing high order Butterworth IIR and EQ filters using the bilinear transform. 

The generated filter coefficients are split out into cascaded biquad sections, for easy use in any biquad or second-order section (SOS) implementation.

### Features
- Lowpass, Highpass, Bandpass and Bandstop IIR & EQ filter design 
- Low and High Shelving filter design
- High order Parametric boost/cut EQ filter design
- Biquad and Biquad Chain implementations (for filtering audio buffers with cascaded biquad sections)
- Compact, readable and well-commented codebase


### Unit tests
As with any good audio signal processing toolkit, there are unit tests that provide basic proof of correctness. There are currently 6 primary test cases that check 113 different conditions.

The unit tests live in `main.cpp` and are written using the excellent and compact [Catch](https://github.com/philsquared/Catch) test framework for C++ & Objective-C.

#### Other filter design repos on GitHub
* Vinnie Falco	[DSP Filters](https://github.com/vinniefalco/DSPFilters)
* Mofr	[Butterworth](https://github.com/mofr/Butterworth)

### Terms and Conditions
Over the course of my work writing audio plugins and music apps, I've written a few different IIR and EQ filter design classes (e.g. RBJ cookbook). These classes are a fresh rewrite and should be considered unoptimized reference code, with an emphasis on clarity (and for pedagogical reasons being able to refer to the textbook, see the math and understand the code). For this reason, in the hope they can be useful, these classes are provided under [GPL v3](http://www.gnu.org/licenses/gpl.html). 

<br/>
Butterworth Filter Design - Copyright © 2013 [iroro orife](http://github.com/ruohoruotsi)

Source code is provided under [GPL v3](http://www.gnu.org/licenses/gpl.html)