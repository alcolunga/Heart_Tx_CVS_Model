# Heart_Tx_CVS_Model

Copyright (C) 2019 A. L. Colunga, K. G. Kim, N. P. Woodall,   T. F. Dardas, M. S. Olufsen, and B. E. Carlson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, and merge the Software subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

The authors be cited in any work resulting from the use of the Software. The associated published article arXiv:1812.11857

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, OR DATA, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

Intro:

This repository contains a MATLAB code for a mathematical model to analyze RHC deidentified data from patient electronic health records. The example code is inspired by a previously published cardiovascular systems-level model (Smith et al., 2004), ignoring the influence of right and left ventricular interaction, as well as omitting the influence of inheritance of blood flowing through the four heart valves.

.m File Description:

Patient233.m: This file contains the RHC data for Patient 233

DriverBasic.m: Calls load_global.m and solves the ODE using dXdT.m and computes the pressures and volumes. This file also plots a six-panel figure comparing computed pressures and volumes to data.

load_global.m: Set nominal parameters, initial conditions, and upper and lower bounds of ODE model

dXdT.m: This function contains the algebraic and differential expressions that describe the cardiovascular system and model

.mat File Description:

results_opt.mat: optimized parameter values for Patient 233
