# iChimera
[![Build Status](https://github.com/lennoxkeeble/iChimera.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lennoxkeeble/iChimera.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/lennoxkeeble/iChimera/blob/main/LICENSE)
![GitHub last commit](https://img.shields.io/github/last-commit/lennoxkeeble/iChimera)

Julia code which computes gravitational waveforms of extreme-mass-ratio inspirals (EMRIs) based on the "Chimera": a local kludge scheme introduced in <a href="https://arxiv.org/abs/1109.0572">Sopuerta & Yunes, 2011</a>. This implementation modifies the treatment of high-order derivatives by computing them analytically, allowing for computational speedup and improved accuracy. Example usage of the code is provided in the examples folder.

## Authors ##

- [Lennox S. Keeble](https://lennoxkeeble.github.io)
- [Alejandro Cardenas-Avendano](https://www.cardenas.sites.wfu.edu)

## MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.