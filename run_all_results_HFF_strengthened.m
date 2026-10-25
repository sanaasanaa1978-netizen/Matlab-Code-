%% run_all_results_HFF_strengthened.m
% Full manuscript reproduction wrapper.
% 1) Runs the primary Galerkin/continuation/Floquet reproduction driver.
% 2) Recomputes the independent grid-PDE validation, including the higher
%    80 x 39 grid, time-step refinement, a second forcing level, and the
%    limited three-dimensional spanwise secondary-stability calculation.
% 3) Verifies the demanding (6,8) shooting/Floquet refinement and the 3-D
%    sign-change/time-step diagnostics.
% 4) Adds Richardson/grid-order analysis, a deterministic uncertainty budget,
%    and a direct eigenvalue/root resolution of the H=10,J=0.5 parent-table discrepancy.
% 5) Exports the independent-grid Figure 09 and Table V, and assembles
%    the complete generated-table file used for manuscript cross-checking.
%
% Each nested script may clear the workspace, so the wrapper deliberately
% recomputes its own path for every call instead of storing a shared variable.

run(fullfile(fileparts(mfilename('fullpath')),'run_all_results_HFF_final.m'));
run(fullfile(fileparts(mfilename('fullpath')),'independent_grid_validation_HFF.m'));
run(fullfile(fileparts(mfilename('fullpath')),'publication_postprocess_HFF.m'));
run(fullfile(fileparts(mfilename('fullpath')),'export_independent_validation_HFF.m'));
run(fullfile(fileparts(mfilename('fullpath')),'verify_high_order_3D_HFF.m'));
fprintf('\nStrengthened HFF reproduction workflow completed.\n');
