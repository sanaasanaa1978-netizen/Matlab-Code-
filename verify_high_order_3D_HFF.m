%% verify_high_order_3D_HFF.m
% Verification gate for the demanding (Mx,Nz)=(6,8) Galerkin refinement
% and the refined three-dimensional spanwise secondary-stability study.
% Run after run_all_results_HFF_final.m and independent_grid_validation_HFF.m.

root = fileparts(mfilename('fullpath'));
dataDir = fullfile(root,'independent_validation');
names = strings(0,1); values = zeros(0,1); criteria = zeros(0,1); passes = false(0,1);

%% High-order Galerkin shooting / Floquet verification
Tsp = readtable(fullfile(root,'shooting_spatial_convergence_matlab.csv'));
Ttr = readtable(fullfile(root,'transverse_convergence_matlab.csv'));
i68 = find(Tsp.Mx==6 & Tsp.Nz==8,1); i45 = find(Tsp.Mx==4 & Tsp.Nz==5,1);
j68 = find(Ttr.Mx==6 & Ttr.Nz==8,1); j45 = find(Ttr.Mx==4 & Ttr.Nz==5,1);

names(end+1,1)="high_order_68_row_present"; values(end+1,1)=double(~isempty(i68) && ~isempty(j68)); criteria(end+1,1)=1; passes(end+1,1)=~isempty(i68) && ~isempty(j68);
if isempty(i68) || isempty(i45) || isempty(j68) || isempty(j45)
    error('Missing (4,5) or (6,8) rows in convergence files.');
end
perRel = abs(Tsp.period(i68)-Tsp.period(i45))/abs(Tsp.period(i45));
nuRel  = abs(Tsp.Nu(i68)-Tsp.Nu(i45))/abs(Tsp.Nu(i45));
muDiff = abs(Ttr.max_transverse_multiplier_abs(j68)-Ttr.max_transverse_multiplier_abs(j45));

names(end+1,1)="high_order_68_shooting_residual"; values(end+1,1)=Tsp.shooting_residual(i68); criteria(end+1,1)=1e-8; passes(end+1,1)=values(end)<criteria(end);
names(end+1,1)="high_order_68_period_relative_difference"; values(end+1,1)=perRel; criteria(end+1,1)=1e-5; passes(end+1,1)=values(end)<criteria(end);
names(end+1,1)="high_order_68_Nu_relative_difference"; values(end+1,1)=nuRel; criteria(end+1,1)=1e-5; passes(end+1,1)=values(end)<criteria(end);
names(end+1,1)="high_order_68_transverse_multiplier_difference"; values(end+1,1)=muDiff; criteria(end+1,1)=1e-5; passes(end+1,1)=values(end)<criteria(end);
names(end+1,1)="high_order_68_transverse_unstable"; values(end+1,1)=Ttr.max_transverse_multiplier_abs(j68); criteria(end+1,1)=1; passes(end+1,1)=values(end)>criteria(end);
names(end+1,1)="high_order_68_translation_neutral_error"; values(end+1,1)=abs(Ttr.translation_neutral_multiplier_abs(j68)-1); criteria(end+1,1)=2e-5; passes(end+1,1)=values(end)<criteria(end);

%% Locate independent-grid 3-D outputs
spanCandidates = {fullfile(dataDir,'limited_3D_spanwise_stability_matlab.csv'), ...
                  fullfile(root,'limited_3D_spanwise_stability_matlab.csv')};
dtCandidates   = {fullfile(dataDir,'limited_3D_timestep_refinement_matlab.csv'), ...
                  fullfile(root,'limited_3D_timestep_refinement_matlab.csv')};
spanFile=''; dtFile='';
for kk=1:numel(spanCandidates)
    if isfile(spanCandidates{kk}), spanFile=spanCandidates{kk}; break; end
end
for kk=1:numel(dtCandidates)
    if isfile(dtCandidates{kk}), dtFile=dtCandidates{kk}; break; end
end
if isempty(spanFile) || isempty(dtFile)
    error(['The 3-D output files are missing. Run independent_grid_validation_HFF.m first, ' ...
           'or run run_all_results_HFF_strengthened.m, which executes the required steps in order.']);
end
T3 = readtable(spanFile);
Tdt = readtable(dtFile);

requiredVars = {'ky_over_a','growth_exponent_48x23','growth_exponent_64x31','growth_exponent_80x39'};
if ~all(ismember(requiredVars,T3.Properties.VariableNames))
    error('The 3-D stability table does not contain the refined 48x23, 64x31 and 80x39 columns.');
end
if ~all(ismember({'ky_over_a','dt','growth_exponent_80x39'},Tdt.Properties.VariableNames))
    error('The 3-D time-step table must contain the 80x39 near-neutral refinement columns.');
end

requiredKy=[0.15 0.19 0.192 0.194 0.1955 0.1956 0.1957 0.196 0.20];
for kk=1:numel(requiredKy)
    if ~any(abs(T3.ky_over_a-requiredKy(kk))<1e-12)
        error('3-D stability file is incomplete: missing ky/a = %.4g.',requiredKy(kk));
    end
end

getSig = @(col,k) col(find(abs(T3.ky_over_a-k)<1e-12,1));
s48_019  = getSig(T3.growth_exponent_48x23,0.19);
s48_0192 = getSig(T3.growth_exponent_48x23,0.192);
s64_0194 = getSig(T3.growth_exponent_64x31,0.194);
s64_0196 = getSig(T3.growth_exponent_64x31,0.196);
s80_015  = getSig(T3.growth_exponent_80x39,0.15);
s80_01956= getSig(T3.growth_exponent_80x39,0.1956);
s80_01957= getSig(T3.growth_exponent_80x39,0.1957);
s80_020  = getSig(T3.growth_exponent_80x39,0.20);

checks = { ...
  "spanwise_48x23_bracket_0p190_0p192", s48_019>0 && s48_0192<0, s48_019*s48_0192, 0; ...
  "spanwise_64x31_bracket_0p194_0p196", s64_0194>0 && s64_0196<0, s64_0194*s64_0196, 0; ...
  "spanwise_80x39_long_wave_growth_0p15", s80_015>0, s80_015, 0; ...
  "spanwise_80x39_bracket_0p1956_0p1957", s80_01956>0 && s80_01957<0, s80_01956*s80_01957, 0; ...
  "spanwise_80x39_decay_0p20", s80_020<0, s80_020, 0};
for kk=1:size(checks,1)
    names(end+1,1)=checks{kk,1}; passes(end+1,1)=checks{kk,2}; values(end+1,1)=checks{kk,3}; criteria(end+1,1)=checks{kk,4};
end

% Finite-value checks ensure that the intended refinement rows were actually run.
for k=[0.15 0.194 0.196 0.20]
    v=getSig(T3.growth_exponent_64x31,k);
    names(end+1,1)="spanwise_64x31_finite_at_"+strrep(sprintf('%.4f',k),'.','p');
    values(end+1,1)=v; criteria(end+1,1)=0; passes(end+1,1)=isfinite(v);
end
for k=[0.15 0.192 0.194 0.1955 0.1956 0.1957 0.196 0.20]
    v=getSig(T3.growth_exponent_80x39,k);
    names(end+1,1)="spanwise_80x39_finite_at_"+strrep(sprintf('%.4f',k),'.','p');
    values(end+1,1)=v; criteria(end+1,1)=0; passes(end+1,1)=isfinite(v);
end

% Near-neutral time-step consistency at ky/a=0.1956 on the finest grid.
mask=abs(Tdt.ky_over_a-0.1956)<1e-12;
v=Tdt.growth_exponent_80x39(mask);
if numel(v)<2, error('Need at least two 80x39 time steps at ky/a=0.1956.'); end
spread=max(v)-min(v);
names(end+1,1)="spanwise_80x39_dt_spread_0p1956"; values(end+1,1)=spread; criteria(end+1,1)=1e-6; passes(end+1,1)=spread<criteria(end);

Tout=table(names(:),values(:),criteria(:),passes(:),'VariableNames',{'check','value','criterion','pass'});
writetable(Tout,fullfile(root,'matlab_extension_verification.csv'));

% Export the demanding (6,8) result from the generated convergence tables.
provenance = "generated by verify_high_order_3D_HFF.m from convergence outputs";
T68=table(Tsp.Mx(i68),Tsp.Nz(i68),0.2,Tsp.period(i68),Tsp.omega(i68),Tsp.r(i68),Tsp.Nu(i68), ...
    Tsp.shooting_residual(i68),Ttr.translation_neutral_multiplier_abs(j68), ...
    Ttr.max_transverse_multiplier_abs(j68),Ttr.transverse_exponent(j68),provenance, ...
    'VariableNames',{'Mx','Nz','deltaRa','period','omega','r','Nu','shooting_residual', ...
    'translation_neutral_multiplier_abs','max_transverse_multiplier_abs','transverse_exponent','provenance'});
writetable(T68,fullfile(root,'high_order_68_verified.csv'));

fid=fopen(fullfile(root,'MATLAB_EXTENSION_VERIFICATION.txt'),'w');
fprintf(fid,'LTNE HFF high-order / refined-3D verification\n');
fprintf(fid,'Generated by verify_high_order_3D_HFF.m\n\n');
for k=1:height(Tout)
    state='FAIL'; if Tout.pass(k), state='PASS'; end
    fprintf(fid,'%-52s : %s   value = %.12g   criterion = %.12g\n', ...
        char(Tout.check(k)),state,Tout.value(k),Tout.criterion(k));
end
overall=all(Tout.pass(:)); state='FAIL'; if overall, state='PASS'; end
fprintf(fid,'\nOVERALL_EXTENSION_PASS = %s\n',state);
fclose(fid);

if ~overall
    error('High-order / refined-3D verification failed. See MATLAB_EXTENSION_VERIFICATION.txt');
end
fprintf('\nHigh-order / refined-3D verification: PASS\n');
