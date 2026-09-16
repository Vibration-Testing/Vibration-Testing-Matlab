% Modal_Analysis
% Version 2026.09.16 16-Sep-2026
% Modal parameter identification, curve fitting, and related utilities.
%
% Curve fitting / mode identification
%   mdofcf     - Curve fit to a multi-degree-of-freedom (MDOF) FRF.
%   sdofcf     - Curve fit to a single-degree-of-freedom (SDOF) FRF.
%   sdof       - SDOF curve-fitting method (Allemang, p. 6-77).
%   mmcf       - MDOF curve fit using the polyreference method.
%   mmcfnr     - MDOF curve fit (non-recursive variant of mmcf).
%   polyref    - MDOF curve fit using polyreference.
%   polyref2   - MDOF curve fit using polyreference (alternate variant).
%   cfitex     - Curve-fit example/demonstration.
%   cferror    - Curve-fit error evaluation.
%   modefit    - Mode-shape fitting utility.
%
% Mode indicators / correlation
%   cmif       - Complex Mode Indicator Function (CMIF), calculates and plots.
%   cmifts     - CMIF test script.
%   mac        - Modal Assurance Criterion between two mode shape sets.
%   comac      - Coordinate Modal Assurance Criterion.
%
% Simulation / model generation
%   frfgen2    - Generates synthetic FRF data.
%   femgui     - GUI for simple finite-element model generation.
%
% Scripts / examples (not reusable functions)
%   makehelp, Uscript, Vscript, sdtest, sec8p1, vtb7_4_1, vtb74test,
%   repmodestest, single_thesis_modes_simple - development, textbook-
%   example, and demonstration scripts.
%
% archive/
%   Superseded/experimental variants kept for reference only
%   (mdofcftry, mdofcfuse, mmcftry, sdofcf2, sdofcf2a, sdofcfold,
%   sdofcvscript); prefer the functions listed above.
