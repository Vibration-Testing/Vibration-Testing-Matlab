% Vibration-Testing-Matlab
% Version 2026.09.16 16-Sep-2026
% Vibration testing, signal-processing, and modal-analysis functions
% affiliated with the manuscript "Vibration Testing with Modal Analysis
% and Health Monitoring - MATLAB version".
%
% Signal processing / spectral estimation
%   asd       - One-sided auto (power) spectral density estimate.
%   crsd      - One-sided cross spectral density estimate (Sxy).
%   coh       - Coherence between two signals.
%   crcor     - Cross correlation.
%   irf       - Impulse response function estimate from a frequency response.
%   tfest     - Estimates a transfer (frequency response) function.
%   tfplot    - Plots magnitude/phase of transfer functions.
%   tfplot2   - Plots transfer functions (alternate style).
%
% Windows
%   blackwin  - Blackman window.
%   boxwin    - Boxcar (rectangular) window.
%   expwin    - N-point exponential window.
%   hammwin   - Hamming window.
%   hannwin   - Von Hann ("hanning") window (calls vonhann).
%   vonhann   - Von Hann window; also applies a window to data.
%   parzen    - Parzen window.
%   triwin    - Bartlett (triangular) window.
%
% Modal / curve fitting
%   sdofcf    - Curve fit to a single-degree-of-freedom (SDOF) FRF.
%   mdofcf    - Curve fit to a multi-degree-of-freedom (MDOF) FRF.
%
% Second-order system utilities
%   so2ss     - Convert second-order (M, C, K) system to state-space form.
%   soeig     - Natural frequencies and mode shapes of a second-order system.
%   sostf     - Transfer function from second-order system matrices.
%   ssim      - Linear second-order system time simulation.
%   ssit      - Subspace iteration eigensolver.
%   guyanold  - Guyan reduction of a second-order system.
%   household - Generates a Householder reflection matrix.
%
% See also the Modal_Analysis/ subfolder and the demos/ subfolder for
% example scripts.
