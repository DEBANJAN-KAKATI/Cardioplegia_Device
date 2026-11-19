%RUN_EVERYTHING  One-click: steady + transient validation and all cases.
%   Results (figures, CSV tables, .mat, logs) are written to ../results
%   Total run time ~20 min on a laptop.
clear; clc;
V  = cp_validate();             % steady verification (analytical solutions)
T  = run_all_cases();           % steady cases 4(a)-(i), design maps
VT = cp_validate_transient();   % transient verification
R  = run_transient_cases();     % digital-twin scenarios D1-D5
