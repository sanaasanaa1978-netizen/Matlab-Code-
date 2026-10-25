%% export_independent_validation_HFF.m
% Reproducible export of the independent-grid manuscript products.
% Requires independent_grid_validation_HFF.m to have completed successfully.
% It generates:
%   Figure_09.{fig,eps,pdf,png,svg}
%   Table_V_gridvalidation_generated.tex (independent-grid/refined-3D table)
%   tables_generated.tex (all six final computational-result tables)

root=fileparts(mfilename('fullpath'));
outdir=fullfile(root,'independent_validation');
transitionFile=fullfile(outdir,'independent_transition_64x31_matlab.csv');
gridFile=fullfile(outdir,'independent_grid_convergence_matlab_with_reference.csv');
spanFile=fullfile(outdir,'limited_3D_spanwise_stability_matlab.csv');
coreTableFile=fullfile(root,'tables_core_generated.tex');
required={transitionFile,gridFile,spanFile,coreTableFile};
for k=1:numel(required)
    if ~isfile(required{k})
        error('Required generated file not found: %s. Run run_all_results_HFF_final.m and independent_grid_validation_HFF.m first.',required{k});
    end
end

%% Figure 09: independent odd-parity growth
D=readtable(transitionFile);
t=D.t; y=D.odd_L2_norm;
omega_c=1.35746754530479; % benchmark value verified by matlab_run_verification.csv
Tcrit=2*pi/omega_c;
fitMask=t<=8*Tcrit;
pp=polyfit(t(fitMask),log(y(fitMask)),1);
sigmaFit=pp(1);

fig=figure('Color','w','Position',[100 100 900 570]);
semilogy(t,y,'LineWidth',1.55,'Color',[0.0000 0.4470 0.6980]); hold on;
tf=linspace(min(t(fitMask)),max(t(fitMask)),250).';
semilogy(tf,exp(polyval(pp,tf)),'--','LineWidth',1.8,'Color',[0.8353 0.3686 0.0000]);
xlabel('Time','Interpreter','latex');
ylabel('Odd-parity $L^2$ norm','Interpreter','latex');
legend({'Independent grid: odd-parity norm',sprintf('Eight-period fit, slope=%.3f',sigmaFit)}, ...
       'Location','northwest','Interpreter','tex','Box','off');
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',0.9,'TickDir','out');
grid on; box on;
saveFiveIndependent(fig,root,'Figure_09',false);
close(fig);

fid=fopen(fullfile(root,'Figure_09_FIT_METADATA.txt'),'w');
fprintf(fid,'Figure 09 source: independent_validation/independent_transition_64x31_matlab.csv\n');
fprintf(fid,'Fit interval: 0 <= t <= 8*(2*pi/omega_c), omega_c = %.15g\n',omega_c);
fprintf(fid,'Fit end time: %.12g\n',8*Tcrit);
fprintf(fid,'Log-linear least-squares slope: %.12g\n',sigmaFit);
fprintf(fid,'Galerkin transverse Floquet exponent: %.12g\n',0.188962);
fprintf(fid,'Relative difference: %.6f%%\n',100*abs(sigmaFit-0.188962)/0.188962);
fclose(fid);

%% Independent-grid computational table: grid + refined limited 3-D stability
G=readtable(gridFile);
S=readtable(spanFile);
fid=fopen(fullfile(root,'Table_V_gridvalidation_generated.tex'),'w');
if fid<0, error('Unable to create Table_V_gridvalidation_generated.tex'); end
fprintf(fid,'\\begin{table}[t]\n\\centering\n');
fprintf(fid,'\\caption{Independent grid refinement at $\\mathrm{Ra}-\\mathrm{Ra}_c=1$ and limited three-dimensional spanwise stability of the resulting traveling state. Positive $\\sigma_{3D}$ denotes growth. \\textit{Source: Author''s calculations.}}\n');
fprintf(fid,'\\label{tab:gridvalidation}\n\\small\n');
fprintf(fid,'\\begin{tabular}{lrrrrrr}\n\\toprule\n');
fprintf(fid,'resolution & $\\Delta t$ & $\\overline{Nu}_{\\rm SW}$ & $|c_{\\rm TW}|$ & $\\omega_{\\rm TW}$ & $\\overline{Nu}_{\\rm TW}$ & $|\\Theta_{11}|$\\\\\n\\midrule\n');
for k=1:height(G)-1
    fprintf(fid,'$%d\\times%d$ & %.4f & %.6f & %.6f & %.6f & %.6f & %.6f\\\\\n', ...
        G.Nx(k),G.Nz(k),G.dt(k),G.Nu_standing(k),G.wave_speed_abs(k),G.omega_traveling(k),G.Nu_traveling(k),G.theta11_amp(k));
end
k=height(G);
fprintf(fid,'Galerkin & -- & %.6f & %.6f & %.6f & %.6f & %.6f\\\\\n', ...
    G.Nu_standing(k),G.wave_speed_abs(k),G.omega_traveling(k),G.Nu_traveling(k),G.theta11_amp(k));
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\n\\vspace{0.6em}\n');
fprintf(fid,'\\begin{tabular}{rrrr}\n\\toprule\n');
fprintf(fid,'$k_y/a$ & $\\sigma_{3D}$, $48\\times23$ & $\\sigma_{3D}$, $64\\times31$ & $\\sigma_{3D}$, $80\\times39$\\\\\n\\midrule\n');
keep=[0.15 0.18 0.192 0.194 0.1955 0.1956 0.1957 0.196 0.20];
for ik=1:numel(keep)
    k=find(abs(S.ky_over_a-keep(ik))<1e-12,1);
    if isempty(k), error('Missing required spanwise row ky/a=%g',keep(ik)); end
    vals=[S.growth_exponent_48x23(k),S.growth_exponent_64x31(k),S.growth_exponent_80x39(k)];
    txt=cell(1,3);
    for j=1:3
        if isnan(vals(j))
            txt{j}='--';
        elseif vals(j)==0
            txt{j}='$0$';
        else
            e=floor(log10(abs(vals(j)))); m=vals(j)/10^e;
            txt{j}=sprintf('$% .4f\\times10^{%d}$',m,e);
        end
    end
    fprintf(fid,'%.4f & %s & %s & %s\\\\\n',S.ky_over_a(k),txt{1},txt{2},txt{3});
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
fclose(fid);

%% Assemble all six final tables at the explicit insertion marker.
core=fileread(coreTableFile);
grid=fileread(fullfile(root,'Table_V_gridvalidation_generated.tex'));
% MATLAB fprintf('%% ...') writes a SINGLE literal percent sign.
% The core producer therefore emits '% __GRID_TABLE_INSERTION_POINT__',
% although older/manually-edited table files may contain '%% ...'.
% Accept either complete marker, preferring the two-percent version to avoid
% accidentally leaving a stray percent in the generated LaTeX.
markerSingle='% __GRID_TABLE_INSERTION_POINT__';
markerDouble='%% __GRID_TABLE_INSERTION_POINT__';
if contains(core,markerDouble)
    marker=markerDouble;
else
    marker=markerSingle;
end
if count(core,marker)~=1
    error(['Expected exactly one grid-table insertion marker in tables_core_generated.tex. ' ...
           'Found %d. Regenerate that file using run_all_results_HFF_final.m, ' ...
           'or inspect the line containing __GRID_TABLE_INSERTION_POINT__.'], ...
           count(core,marker));
end
combined=strrep(core,marker,grid);
outputFile=fullfile(root,'tables_generated.tex');
fid=fopen(outputFile,'w');
if fid<0, error('Unable to create %s',outputFile); end
fwrite(fid,combined);
fclose(fid);
fprintf('Independent manuscript exports complete. Figure 09 fit slope = %.9f\n',sigmaFit);

function saveFiveIndependent(fig,root,name,rasterHeavy)
    if nargin<4, rasterHeavy=false; end
    drawnow;
    % The older MATLAB print(...,'-dpdf') fallback uses a default paper page
    % and can crop a wide figure. Give the PDF a page that fits the figure.
    previousUnits=get(fig,'Units');
    set(fig,'Units','inches');
    figPosition=get(fig,'Position');
    set(fig,'Units',previousUnits);
    margin=0.25;
    figWidth=max(figPosition(3),1);
    figHeight=max(figPosition(4),1);
    set(fig,'PaperUnits','inches', ...
        'PaperSize',[figWidth+2*margin,figHeight+2*margin], ...
        'PaperPosition',[margin,margin,figWidth,figHeight], ...
        'PaperPositionMode','manual');
    savefig(fig,fullfile(root,[name '.fig']));
    if rasterHeavy
        print(fig,fullfile(root,[name '.eps']),'-depsc','-opengl','-r600');
        if exist('exportgraphics','file')==2
            exportgraphics(fig,fullfile(root,[name '.pdf']),'ContentType','image','Resolution',600);
            exportgraphics(fig,fullfile(root,[name '.png']),'Resolution',600,'ContentType','image');
        else
            print(fig,fullfile(root,[name '.pdf']),'-dpdf','-opengl','-r600');
            print(fig,fullfile(root,[name '.png']),'-dpng','-opengl','-r600');
        end
        print(fig,fullfile(root,[name '.svg']),'-dsvg','-opengl');
    else
        print(fig,fullfile(root,[name '.eps']),'-depsc','-painters');
        if exist('exportgraphics','file')==2
            exportgraphics(fig,fullfile(root,[name '.pdf']),'ContentType','vector');
            exportgraphics(fig,fullfile(root,[name '.png']),'Resolution',600);
        else
            print(fig,fullfile(root,[name '.pdf']),'-dpdf','-painters');
            print(fig,fullfile(root,[name '.png']),'-dpng','-painters','-r600');
        end
        print(fig,fullfile(root,[name '.svg']),'-dsvg','-painters');
    end
end
