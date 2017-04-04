# pycs2
Python translations of Cybership II algorithms developed while working on my thesis
[Lindegaard (2003)](https://brage.bibsys.no/xmlui/handle/11250/259448).

## Why Python?
 Matlab was everyone's main tool back in 1999-2001, and essentialy every module
 was drafted and tested in Matlab before porting them all to C  or C++ as Simulink
 S-functions that could be compiled and executed on Cybership II. It worked OK, 
 but the code is hardly readable for anyone, not even me, anymore.
 
 However, I hope that the re-implementation in Python using standard libs like
 Numpy will be more accessible and encourage you to use and entend whatever I 
 have to your projects and needs.

## Installation and requirements

 Python 2.7 and 3.X are supported. Dependencies are pretty basic (so far); [Numpy](http://www.numpy.org/)
 
 Tip: Use *pip* to install dependencies automatically.
 
 ```
 pip install -r requirements.txt
 ```
 
# Algorithms

 This section covers the different parts comprising the final control system on 
 Cybership II. Each module is given its own subfolder
 
 ## rudderthralloc
 
  In [Lindegaard & Fossen (2002)](http://ieeexplore.ieee.org/document/1255661/) we
  developed a novel thrust allocation algorithm for zero/low speed applications
  like dynamic positioning operations. The main objective was to make a continuous
  mapping from the commanded thrust over to propeller setpoints and rudder
  deflection angles.
  
  The mapping consists of two steps.
  
   1. ForceAllocation - map commanded thrust to generalized forces.
   1. Force2Setpoint - map generalized forces to propeller/rudder setpoints.
   
  Currently, only step 1 has been re-implemented.
  
  For examples of use, please see the unit test file *test_forcealloc.py*.
 
# License

 BSD 2-Clause License

