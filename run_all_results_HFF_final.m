%% run_all_results_HFF_final.m
% Complete reproducibility script for:
% Hopf continuation, Floquet stability, O(2)-type symmetry breaking,
% traveling-wave selection, Lyapunov-coefficient maps, steady-Hopf
% interaction, and finite-amplitude energy-gap screening for LTNE Darcy
% convection with Jeffreys thermal memory.
%
% The script saves each figure separately in FIG, EPS, PDF, PNG (600 dpi),
% and SVG formats in the same folder as this file.
%
% No optional toolbox is required.

clear; close all; clc;
root = fileparts(mfilename('fullpath'));
if isempty(root), root = pwd; end

% Demanding one-level-higher Galerkin check requested for final review.
% Set true only when sufficient runtime is available. Results must not be
% reported in the manuscript until this branch has actually completed.
RUN_DEMANDING_HIGH_ORDER = true;

A = 0.13164;
epsPor = 0.30;
gamma = epsPor/(1-epsPor);

C = publicationColors();
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',14);
set(groot,'defaultAxesLineWidth',1.0);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLineLineWidth',2.0);
set(groot,'defaultLineMarkerSize',7);
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLegendFontSize',12);
set(groot,'defaultTextInterpreter','latex');

%% 1. Linear validation against the parent paper
parent = [ ...
    1,  1.00, 49.8737, 41.4139, 10.30; ...
    1,  0.50, 98.5000, 41.4139, 10.30; ...
    1,  0.10,486.3000, 41.4139, 10.30; ...
    5,  1.00,  3.4259, 48.1496, 11.73; ...
    5,  0.50,  5.8302, 48.1496, 11.73; ...
    5,  0.25, 10.6004, 48.1496, 11.73; ...
    5,  0.10, 24.7780, 48.1496, 11.73; ...
   10,  1.00,  1.7499, 55.0441, 12.80; ...
   10,  0.50,  2.5010, 55.0441, 12.80; ...
   10,  0.10,  8.3610, 55.0441, 12.80];

nval = size(parent,1);
computedRatio=zeros(nval,1); computedRstat=zeros(nval,1);
computedAstat2=zeros(nval,1); ratioErrPct=zeros(nval,1);
for k=1:nval
    H=parent(k,1); J=parent(k,2);
    [Mt,Rstat,aStat2]=transitionM(H,J,A,gamma);
    computedRatio(k)=Mt/J;
    computedRstat(k)=Rstat;
    computedAstat2(k)=aStat2;
    ratioErrPct(k)=100*(computedRatio(k)-parent(k,3))/parent(k,3);
end
Tval=table(parent(:,1),parent(:,2),parent(:,3),computedRatio,ratioErrPct, ...
    parent(:,4),computedRstat,parent(:,5),computedAstat2, ...
    'VariableNames',{'H','J','parent_ratio','computed_ratio','ratio_error_pct', ...
    'parent_Rastat','computed_Rastat','parent_astat2','computed_astat2'});
writetable(Tval,fullfile(root,'validation_parent_matlab.csv'));

%% 1b. Exact Fourier/no-memory and interphase-coupling limits
Hlim=10; [Rfour,a2four]=critStat(Hlim,gamma);
MJ=[0.25;0.5;2;10]; limRows=zeros(numel(MJ),6);
for k=1:numel(MJ)
    Lf=LmodalBase(Rfour,a2four,Hlim,MJ(k),MJ(k),1,1,A,gamma);
    ev=sort(real(eig(Lf)),'descend');
    limRows(k,:)=[MJ(k),Rfour,sqrt(a2four),ev(1),ev(2),ev(3)];
end
Tfour=array2table(limRows,'VariableNames',{'M_equals_J','Ra_fourier','a_fourier', ...
    'eig1','eig2','eig3'});
writetable(Tfour,fullfile(root,'fourier_limit_validation_matlab.csv'));
Htest=[1e-6;1e6]; crows=zeros(2,5);
for k=1:2
    [Rtmp,a2tmp]=critStat(Htest(k),gamma);
    if k==1, Rexact=4*pi^2; else, Rexact=4*pi^2*(1+1/gamma); end
    crows(k,:)=[Htest(k),Rtmp,sqrt(a2tmp),Rexact,100*(Rtmp-Rexact)/Rexact];
end
Tcouple=array2table(crows,'VariableNames',{'H','Ra_numeric','a_numeric','Ra_exact_limit','relative_error_pct'});
writetable(Tcouple,fullfile(root,'coupling_limit_validation_matlab.csv'));

%% 2. Smooth linear and normal-form families (Figure 1 plus diagnostics)
families=[1,1,2.0; 5,1,3.0; 10,0.5,4.0];
styles={'-','--','-.'}; marks={'o','s','^'};
F=struct([]);
for jf=1:size(families,1)
    H=families(jf,1); J=families(jf,2); smax=families(jf,3);
    [Mt,Rstat]=transitionM(H,J,A,gamma);
    sn=linspace(1.05,smax,80).';
    Np=numel(sn);
    Raosc=zeros(Np,1); ac2=zeros(Np,1); omega=zeros(Np,1);
    ell1=zeros(Np,1); kappa=zeros(Np,1);
    for k=1:Np
        M=sn(k)*Mt;
        [Raosc(k),ac2(k),omega(k)]=critHopf(H,J,M,A,gamma);
        nf=hopfNormalForm(H,J,M,Raosc(k),ac2(k),A,gamma);
        ell1(k)=nf.ell1; kappa(k)=nf.kappa;
    end
    F(jf).H=H; F(jf).J=J; F(jf).Mt=Mt; F(jf).Rstat=Rstat;
    F(jf).sn=sn; F(jf).Raosc=Raosc; F(jf).ac2=ac2;
    F(jf).omega=omega; F(jf).ell1=ell1; F(jf).kappa=kappa;
end

fig=figure('Position',[100 100 780 470]); hold on; grid on; box on;
for jf=1:numel(F)
    plot(F(jf).sn,F(jf).Raosc/F(jf).Rstat,styles{jf}, ...
        'Color',C(jf,:),'Marker',marks{jf},'MarkerFaceColor','none', ...
        'MarkerIndices',1:10:numel(F(jf).sn), ...
        'DisplayName',sprintf('$H=%g,\\ J=%g$',F(jf).H,F(jf).J));
end
xlabel('Normalized memory ratio $\mathcal{M}$'); ylabel('$Ra_{\rm osc}/Ra_{\rm stat}$');
styleAxes(gca); placeLegend(gca,'northeast'); saveFive(fig,root,'Figure_01',false);

fig=figure('Position',[100 100 780 470]); hold on; grid on; box on;
for jf=1:numel(F)
    plot(F(jf).sn,F(jf).ell1,styles{jf}, ...
        'Color',C(jf,:),'Marker',marks{jf},'MarkerFaceColor','none', ...
        'MarkerIndices',1:10:numel(F(jf).sn), ...
        'DisplayName',sprintf('$H=%g,\\ J=%g$',F(jf).H,F(jf).J));
end
yline(0,'-','Color',[0.35 0.35 0.35],'LineWidth',1.2,'HandleVisibility','off');
xlabel('Normalized memory ratio $\mathcal{M}$'); ylabel('First Lyapunov coefficient $\ell_1$');
styleAxes(gca); placeLegend(gca,'northwest'); saveFive(fig,root,'Diagnostic_Hopf_ell1_families',false);

fig=figure('Position',[100 100 780 470]); hold on; grid on; box on;
for jf=1:numel(F)
    plot(F(jf).sn,F(jf).omega,styles{jf}, ...
        'Color',C(jf,:),'Marker',marks{jf},'MarkerFaceColor','none', ...
        'MarkerIndices',1:10:numel(F(jf).sn), ...
        'DisplayName',sprintf('$H=%g,\\ J=%g$',F(jf).H,F(jf).J));
end
xlabel('Normalized memory ratio $\mathcal{M}$'); ylabel('Critical frequency $\omega_c$');
styleAxes(gca); placeLegend(gca,'northwest'); saveFive(fig,root,'Diagnostic_Hopf_frequency_families',false);

%% 3. Two-dimensional map of ell_1 (Figure 2 plus diagnostics)
Hlist=[1 5 10];
Jgrid=linspace(0.1,1.5,12);
Sgrid=linspace(1.05,4.0,24);
mapRows=[];
for ih=1:numel(Hlist)
    H=Hlist(ih);
    Z=zeros(numel(Sgrid),numel(Jgrid));
    KZ=zeros(size(Z));
    for j=1:numel(Jgrid)
        J=Jgrid(j);
        [Mt,~,~]=transitionM(H,J,A,gamma);
        for is=1:numel(Sgrid)
            M=Sgrid(is)*Mt;
            [Rac,ac2,~]=critHopf(H,J,M,A,gamma);
            nf=hopfNormalForm(H,J,M,Rac,ac2,A,gamma);
            Z(is,j)=nf.ell1; KZ(is,j)=nf.kappa;
            mapRows=[mapRows; H,J,Mt,M,Sgrid(is),Rac,ac2,nf.omega,nf.ell1,nf.kappa]; %#ok<AGROW>
        end
    end
    fig=figure('Position',[100 100 780 490]);
    contourf(Jgrid,Sgrid,Z,24,'LineColor','none'); hold on;
    contour(Jgrid,Sgrid,Z,8,'LineColor',[0.18 0.18 0.18],'LineWidth',0.55);
    cb=colorbar; ylabel(cb,'$\ell_1$','Interpreter','latex','FontSize',14);
    set(cb,'TickLabelInterpreter','latex','FontSize',14);
    xlabel('Jeffreys coefficient $J$'); ylabel('$M/M_t(J,H)$');
    styleAxes(gca); colormap(academicSequentialMap(256));
    caxis([min(Z(:)) max(Z(:))]);
    if H==10
        outName='Figure_02';
    else
        outName=sprintf('Diagnostic_ell1_map_H%g',H);
    end
    saveFive(fig,root,outName,true);
end
Tmap=array2table(mapRows,'VariableNames', ...
    {'H','J','Mt','M','M_over_Mt','Rac','ac2','omega','ell1','kappa'});
writetable(Tmap,fullfile(root,'lyapunov_map_matlab.csv'));

%% 4. Benchmark Hopf point and shooting continuation: Figures 3-6
H=10; J=0.5; M=2;
[Rac0,ac20,omega0]=critHopf(H,J,M,A,gamma);
a0=sqrt(ac20);
nf0=hopfNormalForm(H,J,M,Rac0,ac20,A,gamma);
q=nf0.q; p=nf0.p; G21=nf0.G21; kappa=nf0.kappa; omegaRa=nf0.omegaRa;

deltaRa=[0.02 0.20 0.50 0.80 1.10 1.40 1.70 2.00].';
nB=numel(deltaRa);
rShoot=zeros(nB,1); rWNL=zeros(nB,1); omShoot=zeros(nB,1); omWNL=zeros(nB,1);
NuShoot=zeros(nB,1); phaseShoot=zeros(nB,1); DTShoot=zeros(nB,1); NuHalfShoot=zeros(nB,1);
muSym=zeros(nB,1); shootRes=zeros(nB,1); newtonIts=zeros(nB,1);
xShoot=cell(nB,1); TShoot=zeros(nB,1);

for k=1:nB
    dRa=deltaRa(k);
    model=buildModelGeneric(2,3,a0,H,J,M,Rac0+dRa,A,epsPor,false);
    rWNL(k)=sqrt(-2*kappa*dRa/real(G21));
    omWNL(k)=omega0+omegaRa*dRa+0.5*imag(G21)*rWNL(k)^2;
    X0=zeros(model.N,3); X0(model.idx11c,:)=2*rWNL(k)*real(q).';
    [ts,Ys,~]=integrateSBDF2(model,X0,0.004,80,0.02);
    [xg,Tg]=estimateOrbitGuess(ts,Ys,model,p,20);
    [xs,Ts,Phi,res,it]=shootCorrect(model,xg,Tg,7,1e-10);
    dg=orbitDiagnostics(model,xs,Ts,p);
    muv=eig(Phi); [~,itri]=min(abs(muv-1)); muv(itri)=[];
    rShoot(k)=dg.r; omShoot(k)=2*pi/Ts; NuShoot(k)=dg.Nu;
    phaseShoot(k)=dg.phaseGap; DTShoot(k)=dg.DT; NuHalfShoot(k)=dg.NuHalfAmp;
    muSym(k)=max(abs(muv)); shootRes(k)=res; newtonIts(k)=it;
    xShoot{k}=xs; TShoot(k)=Ts;
end
Tcont=table(deltaRa,Rac0+deltaRa,rShoot,rWNL,omShoot,omWNL,NuShoot,phaseShoot,DTShoot,NuHalfShoot, ...
    muSym,shootRes,newtonIts, ...
    'VariableNames',{'deltaRa','Ra','r_shoot','r_WNL','omega_shoot','omega_WNL', ...
    'Nu_mean','phase_gap_rad','DT_rms','Nu_half_amplitude','max_nontrivial_Floquet', ...
    'shooting_residual','Newton_iterations'});
writetable(Tcont,fullfile(root,'periodic_continuation_matlab.csv'));

fig=figure('Position',[100 100 780 470]); hold on; grid on; box on;
plot(deltaRa,rShoot,'-o','Color',C(1,:),'MarkerFaceColor','none','DisplayName','Periodic shooting');
plot(deltaRa,rWNL,'--','Color',C(2,:),'DisplayName','Weakly nonlinear');
xlabel('$Ra-Ra_c$'); ylabel('Center-manifold amplitude $r$');
styleAxes(gca); placeLegend(gca,'northwest'); saveFive(fig,root,'Figure_03',false);

fig=figure('Position',[100 100 780 470]); hold on; grid on; box on;
plot(deltaRa,omShoot,'-o','Color',C(1,:),'MarkerFaceColor','none','DisplayName','Periodic shooting');
plot(deltaRa,omWNL,'--','Color',C(2,:),'DisplayName','Weakly nonlinear');
xlabel('$Ra-Ra_c$'); ylabel('Oscillation frequency');
styleAxes(gca); placeLegend(gca,'southwest'); saveFive(fig,root,'Figure_04',false);

fig=figure('Position',[100 100 780 470]); grid on; box on;
plot(deltaRa,NuShoot-1,'-o','Color',C(1,:),'MarkerFaceColor','none');
xlabel('$Ra-Ra_c$'); ylabel('$\overline{Nu}-1$');
styleAxes(gca); saveFive(fig,root,'Figure_05',false);

fig=figure('Position',[100 100 780 470]); hold on; grid on; box on;
plot(deltaRa,muSym,'-o','Color',C(4,:),'MarkerFaceColor','none');
yline(1,'--','Color',C(3,:),'LineWidth',1.4);
xlabel('$Ra-Ra_c$'); ylabel('Largest nontrivial Floquet multiplier $|\mu_F|$');
styleAxes(gca); saveFive(fig,root,'Figure_06',false);

%% 5. Spatial convergence of shooting solution at deltaRa=0.2
pairs=[2 3;3 4;4 5];
if RUN_DEMANDING_HIGH_ORDER, pairs=[pairs;6 8]; end
spRows=zeros(size(pairs,1),7);
for k=1:size(pairs,1)
    Mx=pairs(k,1); Nz=pairs(k,2); dRa=0.2;
    model=buildModelGeneric(Mx,Nz,a0,H,J,M,Rac0+dRa,A,epsPor,false);
    rw=sqrt(-2*kappa*dRa/real(G21));
    X0=zeros(model.N,3); X0(model.idx11c,:)=2*rw*real(q).';
    [ts,Ys,~]=integrateSBDF2(model,X0,0.0035,80,0.02);
    [xg,Tg]=estimateOrbitGuess(ts,Ys,model,p,20);
    [xs,Ts,~,res,~]=shootCorrect(model,xg,Tg,7,1e-10);
    dg=orbitDiagnostics(model,xs,Ts,p);
    spRows(k,:)=[Mx,Nz,Ts,2*pi/Ts,dg.r,dg.Nu,res];
end
ref=spRows(end,:);
rRel=100*abs(spRows(:,5)-ref(5))/abs(ref(5));
omRel=100*abs(spRows(:,4)-ref(4))/abs(ref(4));
NuRel=100*abs(spRows(:,6)-ref(6))/abs(ref(6));
Tsp=table(spRows(:,1),spRows(:,2),spRows(:,3),spRows(:,4),spRows(:,5),spRows(:,6), ...
    rRel,omRel,NuRel,spRows(:,7), ...
    'VariableNames',{'Mx','Nz','period','omega','r','Nu','r_rel_pct','omega_rel_pct','Nu_rel_pct','shooting_residual'});
writetable(Tsp,fullfile(root,'shooting_spatial_convergence_matlab.csv'));

%% 5b. Fixed-step temporal convergence of the nonlinear periodic state
dconv=0.2; rconv=sqrt(-2*kappa*dconv/real(G21));
dtList=[0.02 0.01 0.005 0.0025].';
tc=zeros(numel(dtList),8);
for k=1:numel(dtList)
    dt=dtList(k);
    model=buildModelGeneric(5,6,a0,H,J,M,Rac0+dconv,A,epsPor,false);
    X0=zeros(model.N,3); X0(model.idx11c,:)=2*rconv*real(q).';
    saveDt=max(0.02,5*dt);
    [ts,Ys,Nus]=integrateSBDF2(model,X0,dt,180,saveDt);
    dg=diagnosticsPeriodic(ts,Ys,Nus,model,p,100);
    tc(k,:)=[dt,dg.r,dg.omega,dg.Nu,dg.phaseGap,dg.DT,dg.NuHalfAmp,numel(ts)];
end
ref=tc(end,:);
rRel=100*abs(tc(:,2)-ref(2))/abs(ref(2));
omRel=100*abs(tc(:,3)-ref(3))/abs(ref(3));
NuRel=100*abs(tc(:,4)-ref(4))/abs(ref(4));
Ttime=table(tc(:,1),tc(:,2),rRel,tc(:,3),omRel,tc(:,4),NuRel,tc(:,5),tc(:,6),tc(:,7),tc(:,8), ...
    'VariableNames',{'dt','r','r_rel_pct','omega','omega_rel_pct','Nu','Nu_rel_pct', ...
    'phase_gap_rad','DT_rms','Nu_half_amplitude','saved_samples'});
writetable(Ttime,fullfile(root,'time_step_convergence_matlab.csv'));

%% 6. Transverse Floquet analysis and traveling-wave selection: Figures 7-9 plus diagnostics
dTW=[0.2 0.5 1.0 2.0].';
Tend=[520 350 300 260].';
nTW=numel(dTW);
muTrans=zeros(nTW,1); muTranslate=zeros(nTW,1); expTrans=zeros(nTW,1); cTW=zeros(nTW,1); NuTW=zeros(nTW,1);
NuSW=zeros(nTW,1); antiGrowth=zeros(nTW,1);
tsBreak=[]; antiBreak=[]; fieldState=[]; fieldModel=[];

for k=1:nTW
    dRa=dTW(k);
    modelC=buildModelGeneric(2,3,a0,H,J,M,Rac0+dRa,A,epsPor,false);
    rw=sqrt(-2*kappa*dRa/real(G21));
    X0=zeros(modelC.N,3); X0(modelC.idx11c,:)=2*rw*real(q).';
    [ts,Ys,~]=integrateSBDF2(modelC,X0,0.004,80,0.02);
    [xg,Tg]=estimateOrbitGuess(ts,Ys,modelC,p,20);
    [xStand,Tstand,~,~,~]=shootCorrect(modelC,xg,Tg,7,1e-10);
    dgs=orbitDiagnostics(modelC,xStand,Tstand,p);
    NuSW(k)=dgs.Nu;

    modelO=buildModelGeneric(2,3,a0,H,J,M,Rac0+dRa,A,epsPor,true);
    yFull=embedCosState(modelC,modelO,xStand);

    [muT,exT,mu0]=transverseFloquet(modelO,yFull,Tstand);
    muTrans(k)=muT; expTrans(k)=exT; muTranslate(k)=mu0;

    Xfull=reshape(yFull,3,modelO.N).';
    Xfull(modelO.idx11s,:)=Xfull(modelO.idx11s,:)+1e-4*real(q).';
    [tt,YY,NuSeries]=integrateSBDF2(modelO,Xfull,0.004,Tend(k),0.04);
    dgTW=travelingDiagnostics(tt,YY,NuSeries,modelO,Tend(k)-80);
    cTW(k)=abs(dgTW.speed); NuTW(k)=dgTW.Nu; antiGrowth(k)=dgTW.growth;

    if abs(dRa-0.2)<1e-12
        tsBreak=tt; antiBreak=dgTW.anti;
    end
    if abs(dRa-1.0)<1e-12
        fieldState=squeeze(YY(end,:,:));
        fieldModel=modelO;
    end
end

Ttrans=table(dTW,muTranslate,muTrans,expTrans,'VariableNames', ...
    {'deltaRa','translation_neutral_multiplier_abs','max_transverse_multiplier_abs','transverse_exponent'});
writetable(Ttrans,fullfile(root,'transverse_floquet_matlab.csv'));
Ttw=table(dTW,antiGrowth,cTW,NuTW,NuSW,'VariableNames', ...
    {'deltaRa','antisymmetric_growth','wave_speed_abs','Nu_traveling','Nu_standing'});
writetable(Ttw,fullfile(root,'traveling_wave_branch_matlab.csv'));

%% 6b. Targeted nonlinear parameter sensitivity at deltaRa = 1
% This small set is intentionally not a full nonlinear survey. It tests
% whether the 2-D standing-to-traveling selection and heat-transfer ordering
% persist when H and the absolute Jeffreys parameters are changed.
sensCases=[10 0.5 2.0; 5 0.5 4.0; 10 1.0 3.0];
sensRows=zeros(size(sensCases,1),9);
for is=1:size(sensCases,1)
    Hs=sensCases(is,1); Js=sensCases(is,2); Ms=sensCases(is,3);
    [Racs,ac2s,~]=critHopf(Hs,Js,Ms,A,gamma); as=sqrt(ac2s);
    nfs=hopfNormalForm(Hs,Js,Ms,Racs,ac2s,A,gamma);
    [Mts,~,~]=transitionM(Hs,Js,A,gamma);
    dRa=1.0;
    modelC=buildModelGeneric(2,3,as,Hs,Js,Ms,Racs+dRa,A,epsPor,false);
    rw=sqrt(-2*nfs.kappa*dRa/real(nfs.G21));
    X0=zeros(modelC.N,3); X0(modelC.idx11c,:)=2*rw*real(nfs.q).';
    [ts,Ys,~]=integrateSBDF2(modelC,X0,0.004,100,0.02);
    [xg,Tg]=estimateOrbitGuess(ts,Ys,modelC,nfs.p,25);
    [xStand,Tstand,~,~,~]=shootCorrect(modelC,xg,Tg,7,1e-10);
    dgs=orbitDiagnostics(modelC,xStand,Tstand,nfs.p);

    modelO=buildModelGeneric(2,3,as,Hs,Js,Ms,Racs+dRa,A,epsPor,true);
    yFull=embedCosState(modelC,modelO,xStand);
    [~,exT,~]=transverseFloquet(modelO,yFull,Tstand);
    Xfull=reshape(yFull,3,modelO.N).';
    Xfull(modelO.idx11s,:)=Xfull(modelO.idx11s,:)+1e-4*real(nfs.q).';
    [tt,YY,NuSeries]=integrateSBDF2(modelO,Xfull,0.004,320,0.04);
    dgTW=travelingDiagnostics(tt,YY,NuSeries,modelO,240);
    sensRows(is,:)=[Hs,Js,Ms,Ms/Mts,Racs,exT,abs(dgTW.speed),dgs.Nu,dgTW.Nu];
end
Tsens=array2table(sensRows,'VariableNames',{'H','J','M','M_over_Mt','Rac', ...
    'transverse_exponent','wave_speed_abs','Nu_standing','Nu_traveling'});
writetable(Tsens,fullfile(root,'nonlinear_parameter_sensitivity_matlab.csv'));

fig=figure('Position',[100 100 780 470]); hold on; grid on; box on;
plot(dTW,muTrans,'-o','Color',C(3,:),'MarkerFaceColor','none');
yline(1,'--','Color',C(1,:),'LineWidth',1.4);
xlabel('$Ra-Ra_c$'); ylabel('Dominant transverse Floquet multiplier $|\mu_\perp|$');
styleAxes(gca); saveFive(fig,root,'Figure_07',false);

fig=figure('Position',[100 100 780 470]); grid on; box on;
semilogy(tsBreak,max(antiBreak,1e-12),'-','Color',C(4,:));
xlabel('Time'); ylabel('Antisymmetric modal norm');
styleAxes(gca); saveFive(fig,root,'Diagnostic_symmetry_breaking_growth',false);

fig=figure('Position',[100 100 780 470]); grid on; box on;
plot(dTW,cTW,'-o','Color',C(3,:),'MarkerFaceColor','none');
xlabel('$Ra-Ra_c$'); ylabel('Traveling-wave speed $|c_{\rm TW}|$');
styleAxes(gca); saveFive(fig,root,'Diagnostic_traveling_wave_speed',false);

fig=figure('Position',[100 100 780 470]); hold on; grid on; box on;
plot(dTW,NuSW-1,'-o','Color',C(2,:),'MarkerFaceColor','none','DisplayName','Standing wave');
plot(dTW,NuTW-1,'--s','Color',C(1,:),'MarkerFaceColor','none','DisplayName','Traveling wave');
xlabel('$Ra-Ra_c$'); ylabel('$\overline{Nu}-1$');
styleAxes(gca); placeLegend(gca,'northwest'); saveFive(fig,root,'Figure_08',false);

[xp,zp,thetaTW,phiTW]=reconstructFields(fieldModel,fieldState,241,161);
fig=figure('Position',[100 100 720 470]);
contourf(xp,zp,thetaTW,31,'LineColor','none'); hold on;
contour(xp,zp,thetaTW,[0 0],'LineColor',[0.20 0.20 0.20],'LineWidth',0.80);
colormap(academicDivergingMap(256));
lim=max(abs(thetaTW(:))); if lim>0, caxis([-lim lim]); end
cb=colorbar; ylabel(cb,'$\theta$','Interpreter','latex','FontSize',14);
set(cb,'TickLabelInterpreter','latex','FontSize',14);
xlabel('$x$'); ylabel('$z$'); axis tight; styleAxes(gca); saveFive(fig,root,'Diagnostic_traveling_temperature_field',true);

fig=figure('Position',[100 100 720 470]);
contourf(xp,zp,phiTW,31,'LineColor','none'); hold on;
contour(xp,zp,phiTW,[0 0],'LineColor',[0.20 0.20 0.20],'LineWidth',0.80);
colormap(academicDivergingMap(256));
lim=max(abs(phiTW(:))); if lim>0, caxis([-lim lim]); end
cb=colorbar; ylabel(cb,'$\phi$','Interpreter','latex','FontSize',14);
set(cb,'TickLabelInterpreter','latex','FontSize',14);
xlabel('$x$'); ylabel('$z$'); axis tight; styleAxes(gca); saveFive(fig,root,'Diagnostic_traveling_solid_temperature',true);

%% 7. Transverse-Floquet convergence at deltaRa=0.2
pairs=[2 3;3 4;4 5];
if RUN_DEMANDING_HIGH_ORDER, pairs=[pairs;6 8]; end
trRows=zeros(size(pairs,1),6);
for k=1:size(pairs,1)
    Mx=pairs(k,1); Nz=pairs(k,2); dRa=0.2;
    modelC=buildModelGeneric(Mx,Nz,a0,H,J,M,Rac0+dRa,A,epsPor,false);
    rw=sqrt(-2*kappa*dRa/real(G21));
    X0=zeros(modelC.N,3); X0(modelC.idx11c,:)=2*rw*real(q).';
    [ts,Ys,~]=integrateSBDF2(modelC,X0,0.0035,80,0.02);
    [xg,Tg]=estimateOrbitGuess(ts,Ys,modelC,p,20);
    [xs,Ts,~,~,~]=shootCorrect(modelC,xg,Tg,7,1e-10);
    modelO=buildModelGeneric(Mx,Nz,a0,H,J,M,Rac0+dRa,A,epsPor,true);
    yF=embedCosState(modelC,modelO,xs);
    [muT,exT,mu0]=transverseFloquet(modelO,yF,Ts);
    trRows(k,:)=[Mx,Nz,Ts,mu0,muT,exT];
end
ref=trRows(end,:);
muRel=100*abs(trRows(:,5)-ref(5))/abs(ref(5));
exRel=100*abs(trRows(:,6)-ref(6))/abs(ref(6));
Ttrc=table(trRows(:,1),trRows(:,2),trRows(:,3),trRows(:,4),trRows(:,5),trRows(:,6),muRel,exRel, ...
    'VariableNames',{'Mx','Nz','period','translation_neutral_multiplier_abs', ...
    'max_transverse_multiplier_abs','transverse_exponent','mu_rel_pct','exponent_rel_pct'});
writetable(Ttrc,fullfile(root,'transverse_convergence_matlab.csv'));

%% 8. Steady-Hopf interaction: Figure 10
inter=steadyHopfInteraction(10,0.5,A,gamma);
Tint=struct2table(inter);
writetable(Tint,fullfile(root,'steady_hopf_interaction_matlab.csv'));

dR=linspace(0,0.6,301).';
bExist=inter.hopf_exist_slope;
bLow=inter.lower_slope; bHigh=inter.upper_slope;
fig=figure('Position',[100 100 780 500]); hold on; box on; grid on;
fill([dR;flipud(dR)],[(-0.18+0*dR);flipud(bExist*dR)],[0.88 0.92 1.00], ...
    'EdgeColor','none','HandleVisibility','off');
fill([dR;flipud(dR)],[(bExist*dR);flipud(bLow*dR)],[0.78 0.85 0.97], ...
    'EdgeColor','none','HandleVisibility','off');
fill([dR;flipud(dR)],[(bLow*dR);flipud(bHigh*dR)],[0.91 0.84 0.95], ...
    'EdgeColor','none','HandleVisibility','off');
fill([dR;flipud(dR)],[(bHigh*dR);flipud(0.18+0*dR)],[1.00 0.88 0.88], ...
    'EdgeColor','none','HandleVisibility','off');
plot(dR,bExist*dR,':','Color',C(3,:),'LineWidth',1.5,'DisplayName','Hopf existence');
plot(dR,bLow*dR,'--','Color',C(2,:),'LineWidth',2,'DisplayName','Hopf transverse stability');
plot(dR,bHigh*dR,'-','Color',C(1,:),'LineWidth',2,'DisplayName','Steady transverse stability');
text(0.39,-0.14,'steady only','Interpreter','latex','FontSize',12);
text(0.36,0.01,'coexistence','Interpreter','latex','FontSize',12);
text(0.39,0.14,'Hopf only','Interpreter','latex','FontSize',12);
xlim([0 0.6]); ylim([-0.18 0.18]);
xlabel('$\delta_R=Ra-Ra_*$'); ylabel('$\delta_M=M-M_*$');
styleAxes(gca); placeLegend(gca,'northwest'); saveFive(fig,root,'Figure_10',false);

%% 9. Finite-amplitude energy-gap screening (diagnostic figure)
[RaE,aE2]=energyThreshold(H,J,M,gamma);
[RaStat,aS2]=critStat(H,gamma);
RaGap=[50.25 51.00 51.75 52.25];
ampList=[0.25 0.5 1 2 4];
atypeName={'energy','hopf','stationary'};
atypeA2=[aE2 ac20 aS2];
gapRows={};
maxFinal=zeros(numel(RaGap),1);
for ir=1:numel(RaGap)
    RaNow=RaGap(ir);
    for ia=1:3
        a2now=atypeA2(ia);
        model=buildModelGeneric(1,2,sqrt(a2now),H,J,M,RaNow,A,epsPor,false);
        L11=LmodalBase(RaNow,a2now,H,J,M,1,1,A,gamma);
        [V,D]=eig(complex(L11)); ev=diag(D);
        [~,ie]=max(real(ev)); q0=V(:,ie); q0=q0/norm(q0);
        q0=q0*exp(-1i*angle(q0(1))); qr=real(q0); qr=qr/norm(qr);
        for ja=1:numel(ampList)
            X0=zeros(model.N,3); X0(model.idx11c,:)=ampList(ja)*qr.';
            [tg,Yg,~]=integrateSBDF2(model,X0,0.01,120,0.1);
            flat=reshape(Yg,size(Yg,1),[]);
            norms=sqrt(sum(flat.^2,2));
            late=norms(tg>=90);
            lm=max(late); fn=norms(end);
            gapRows(end+1,:)={RaNow,atypeName{ia},a2now,ampList(ja),lm,fn}; %#ok<AGROW>
            maxFinal(ir)=max(maxFinal(ir),fn);
        end
    end
end
Tgap=cell2table(gapRows,'VariableNames',{'Ra','wavenumber_type','a2','initial_amplitude','late_max_norm','final_norm'});
writetable(Tgap,fullfile(root,'energy_gap_search_matlab.csv'));

% Higher-resolution extended checks at the most persistent upper-gap state.
extRows={};
for ia=1:2
    a2now=atypeA2(ia);
    model=buildModelGeneric(2,3,sqrt(a2now),H,J,M,52.25,A,epsPor,false);
    L11=LmodalBase(52.25,a2now,H,J,M,1,1,A,gamma);
    [V,D]=eig(complex(L11)); ev=diag(D); [~,ie]=max(real(ev));
    q0=V(:,ie); q0=q0/norm(q0); q0=q0*exp(-1i*angle(q0(1)));
    qr=real(q0); qr=qr/norm(qr);
    X0=zeros(model.N,3); X0(model.idx11c,:)=4*qr.';
    [tg,Yg,~]=integrateSBDF2(model,X0,0.005,300,0.1);
    flat=reshape(Yg,size(Yg,1),[]);
    norms=sqrt(sum(flat.^2,2));
    extRows(end+1,:)={52.25,atypeName{ia},4,2,3,0.005,300,max(norms(tg>=250)),norms(end)}; %#ok<AGROW>
end
Text=cell2table(extRows,'VariableNames',{'Ra','wavenumber_type','initial_amplitude','Mx','Nz','dt','T','late_max_norm','final_norm'});
writetable(Text,fullfile(root,'energy_gap_extended_matlab.csv'));

fig=figure('Position',[100 100 780 470]); hold on; box on; grid on;
semilogy(RaGap,maxFinal,'-o','Color',C(3,:),'MarkerFaceColor','none','DisplayName','Maximum final norm, $t=120$');
semilogy(52.25,max(Text.final_norm),'s','Color',C(1,:),'MarkerFaceColor','none','MarkerSize',8, ...
    'DisplayName','Extended check, $t=300$');
xlabel('$Ra$'); ylabel('Maximum terminal modal norm');
styleAxes(gca); placeLegend(gca,'northwest'); saveFive(fig,root,'Diagnostic_energy_gap_screening',false);

%% 9b. Export the six generated computational-result tables in LaTeX form
writeFinalLatexTables(root,Tval,Tmap,Tcont,Ttrans,Ttw,Tsens,Tsp,Ttime);

%% 10. Automatic numerical self-check and run-verification report
checks={}; pass=[]; value=[]; tolerance=[];
% One-time algebraic check of the vectorized quadratic contraction.
modelCheck=buildModelGeneric(2,3,a0,H,J,M,Rac0+0.2,A,epsPor,false);
yCheck=0.01*sin((1:(3*modelCheck.N)).');
quadErr=norm(nonlinearVector(modelCheck,yCheck)-nonlinearVectorDirect(modelCheck,yCheck),inf);
checks{end+1,1}='vectorized_quadratic_contraction'; value(end+1,1)=quadErr; tolerance(end+1,1)=1e-12; pass(end+1,1)=quadErr<tolerance(end);
checks{end+1,1}='benchmark_Rac'; value(end+1,1)=Rac0; tolerance(end+1,1)=2e-7; pass(end+1,1)=abs(Rac0-52.4889464007)<tolerance(end);
checks{end+1,1}='benchmark_ac2'; value(end+1,1)=ac20; tolerance(end+1,1)=2e-6; pass(end+1,1)=abs(ac20-11.5155206)<tolerance(end);
checks{end+1,1}='benchmark_omega'; value(end+1,1)=omega0; tolerance(end+1,1)=2e-7; pass(end+1,1)=abs(omega0-1.35746754)<tolerance(end);
checks{end+1,1}='benchmark_ell1'; value(end+1,1)=nf0.ell1; tolerance(end+1,1)=2e-6; pass(end+1,1)=abs(nf0.ell1+1.37352050)<tolerance(end);
checks{end+1,1}='all_Hopf_map_ell1_negative'; value(end+1,1)=max(Tmap.ell1); tolerance(end+1,1)=0; pass(end+1,1)=all(Tmap.ell1<0);
checks{end+1,1}='all_Hopf_map_kappa_positive'; value(end+1,1)=min(Tmap.kappa); tolerance(end+1,1)=0; pass(end+1,1)=all(Tmap.kappa>0);
checks{end+1,1}='shooting_residual_max'; value(end+1,1)=max(shootRes); tolerance(end+1,1)=1e-8; pass(end+1,1)=max(shootRes)<tolerance(end);
checks{end+1,1}='symmetric_Floquet_max'; value(end+1,1)=max(muSym); tolerance(end+1,1)=1; pass(end+1,1)=all(muSym<1);
checks{end+1,1}='translation_neutral_error_max'; value(end+1,1)=max(abs(muTranslate-1)); tolerance(end+1,1)=2e-5; pass(end+1,1)=value(end)<tolerance(end);
checks{end+1,1}='transverse_unstable_min'; value(end+1,1)=min(muTrans); tolerance(end+1,1)=1; pass(end+1,1)=all(muTrans>1);
linearPhaseSpeed=omega0/a0;
checks{end+1,1}='near_onset_wave_speed_consistency'; value(end+1,1)=abs(cTW(1)-linearPhaseSpeed); tolerance(end+1,1)=1e-4; pass(end+1,1)=value(end)<tolerance(end);
checks{end+1,1}='energy_gap_extended_final_norm_max'; value(end+1,1)=max(Text.final_norm); tolerance(end+1,1)=1e-6; pass(end+1,1)=value(end)<tolerance(end);
Tverify=table(checks,logical(pass),value,tolerance,'VariableNames',{'check','pass','value','tolerance'});
writetable(Tverify,fullfile(root,'matlab_run_verification.csv'));
fid=fopen(fullfile(root,'MATLAB_RUN_VERIFICATION.txt'),'w');
fprintf(fid,'LTNE Jeffreys advanced calculation verification report\n');
fprintf(fid,'Generated automatically by run_all_results_HFF_final.m\n\n');
for k=1:height(Tverify)
    if Tverify.pass(k), status='PASS'; else, status='FAIL'; end
    fprintf(fid,'%-40s : %s   value = %.12g   criterion = %.12g\n', ...
        Tverify.check{k},status,Tverify.value(k),Tverify.tolerance(k));
end
if all(Tverify.pass), overall='PASS'; else, overall='FAIL'; end
fprintf(fid,'\nOVERALL_PASS = %s\n',overall);
fclose(fid);
if ~all(Tverify.pass)
    error('One or more numerical self-checks failed. Inspect MATLAB_RUN_VERIFICATION.txt.');
end

fprintf('\nAll calculations completed and automatic verification checks passed.\n');
fprintf('Benchmark: Rac = %.10f, ac^2 = %.10f, omega = %.10f, ell1 = %.10f\n', ...
    Rac0,ac20,omega0,nf0.ell1);
fprintf('Energy threshold: RaE = %.10f, aE^2 = %.10f\n',RaE,aE2);
fprintf('Output folder: %s\n',root);

%% ========================================================================
% Local functions
% ========================================================================

function L=LmodalBase(Ra,aBase2,H,J,M,m,n,A,gamma)
    a=sqrt(aBase2);
    Lam=(m*a)^2+(n*pi)^2;
    alpha=(1/M)*(1-J/M);
    if m==0, buoy=0; else, buoy=Ra*(m*a)^2/Lam; end
    L=[buoy-Lam-H, H, 0; ...
       gamma*H/A, -((J/M)*Lam+gamma*H)/A, -alpha*Lam/A; ...
       0, 1, -1/M];
end

function [Rstat,a2]=critStat(H,gamma)
    f=@(x) ((pi^2+x).^2./x)+H*(pi^2+x).^2./(x.*(pi^2+x+gamma*H));
    [a2,Rstat]=fminbnd(f,2,30,optimset('TolX',1e-11,'Display','off'));
end

function Ra=hopfRaAtA2(a2,H,J,M,A,gamma)
    Lam=a2+pi^2;
    alpha=(1/M)*(1-J/M);
    d=-((J/M)*Lam+gamma*H)/A;
    e=-alpha*Lam/A;
    f=-1/M;
    b=gamma*H/A;
    u=d+f;
    w=d*f-e;
    v=w-H*b;
    rr=roots([-u, H*b-u^2, -(u*v+H*b*f)]);
    cand=[];
    for k=1:numel(rr)
        if abs(imag(rr(k)))<1e-9
            x=real(rr(k));
            Rt=(x+Lam+H)*Lam/a2;
            c2=-(x+d+f);
            c1=x*(d+f)+d*f-e-H*b;
            if Rt>0 && c2>0 && c1>0
                cand(end+1)=Rt; %#ok<AGROW>
            end
        end
    end
    if isempty(cand), Ra=NaN; else, Ra=min(cand); end
end

function [Ra,a2,omega]=critHopf(H,J,M,A,gamma)
    ag=linspace(2,30,70); rg=nan(size(ag));
    for k=1:numel(ag), rg(k)=hopfRaAtA2(ag(k),H,J,M,A,gamma); end
    ok=find(isfinite(rg));
    if isempty(ok), error('No admissible Hopf root.'); end
    [~,jj]=min(rg(ok)); x0=ag(ok(jj));
    lo=max(0.2,x0-2); hi=x0+2;
    obj=@(x) finiteHopfObjective(x,H,J,M,A,gamma);
    [a2,Ra]=fminbnd(obj,lo,hi,optimset('TolX',1e-10,'Display','off'));
    L=LmodalBase(Ra,a2,H,J,M,1,1,A,gamma);
    c=poly(L); omega=sqrt(c(3));
end

function y=finiteHopfObjective(x,H,J,M,A,gamma)
    y=hopfRaAtA2(x,H,J,M,A,gamma);
    if ~isfinite(y), y=1e12; end
end

function [Mt,Rstat,aStat2]=transitionM(H,J,A,gamma)
    [Rstat,aStat2]=critStat(H,gamma);
    Ms=logspace(log10(J*1.001),log10(max(200,200*J)),90);
    g=nan(size(Ms));
    for k=1:numel(Ms)
        try
            [Ro,~,~]=critHopf(H,J,Ms(k),A,gamma);
            g(k)=Ro-Rstat;
        catch
            g(k)=NaN;
        end
    end
    Mt=NaN;
    for k=1:numel(Ms)-1
        if isfinite(g(k)) && isfinite(g(k+1)) && g(k)*g(k+1)<0
            Mt=fzero(@(mm) transitionResidual(mm,H,J,A,gamma,Rstat), ...
                [Ms(k),Ms(k+1)],optimset('TolX',1e-9,'Display','off'));
            break;
        end
    end
    if ~isfinite(Mt), error('Transition M could not be bracketed.'); end
end

function g=transitionResidual(M,H,J,A,gamma,Rstat)
    [Ro,~,~]=critHopf(H,J,M,A,gamma);
    g=Ro-Rstat;
end

function nf=hopfNormalForm(H,J,M,Ra,a2,A,gamma)
    L=LmodalBase(Ra,a2,H,J,M,1,1,A,gamma);
    [V,D]=eig(complex(L)); ev=diag(D);
    pos=find(imag(ev)>0); [~,jj]=min(abs(real(ev(pos)))); iq=pos(jj);
    q=V(:,iq); q=q/norm(q); q=q*exp(-1i*angle(q(1)));
    omega=imag(ev(iq));
    [Vp,Dp]=eig(L'); ep=diag(Dp); [~,ip]=min(abs(ep+1i*omega)); p=Vp(:,ip);
    p=p/conj(p'*q);
    Lam=a2+pi^2; chi=pi*sqrt(Ra)*a2/Lam;
    L0=LmodalBase(Ra,a2,H,J,M,0,2,A,gamma);
    bqq=[-chi*q(1)^2;0;0];
    bqb=[-chi*q(1)*conj(q(1));0;0];
    h20=(2i*omega*eye(3)-L0)\bqq;
    h11=L0\bqb;
    fb=[chi*(conj(q(1))*h20(1)-2*q(1)*h11(1));0;0];
    G21=p'*fb;
    dL=zeros(3); dL(1,1)=a2/Lam;
    dlambda=p'*(dL*q);
    nf=struct('q',q,'p',p,'G21',G21,'ell1',real(G21)/(2*omega), ...
        'kappa',real(dlambda),'omegaRa',imag(dlambda),'omega',omega);
end

function model=buildModelGeneric(Mx,Nz,a,H,J,M,Ra,A,epsPor,fullO2)
    gamma=epsPor/(1-epsPor); R=sqrt(Ra);
    modes=[];
    for m=0:Mx
        for n=1:Nz
            modes=[modes; m n 0]; %#ok<AGROW>
            if fullO2 && m>0, modes=[modes; m n 1]; end %#ok<AGROW>
        end
    end
    N=size(modes,1); Lx=2*pi/a;
    Nxq=max(32,4*(2*Mx+1)); Nzq=max(24,3*Nz+8);
    x=(0:Nxq-1)'*Lx/Nxq;
    [z,wz]=gaussLegendre(Nzq,0,1);
    [Xg,Zg]=ndgrid(x,z);
    W2=(Lx/Nxq)*ones(Nxq,1)*wz.'; W2=W2(:).';
    P=numel(Xg);
    E=zeros(N,P); Ex=zeros(N,P); Ez=zeros(N,P);
    U=zeros(N,P); Wv=zeros(N,P);
    Lam=zeros(N,1); norms=zeros(N,1); idx11c=[]; idx11s=[];
    Lfull=zeros(3*N);
    for k=1:N
        m=modes(k,1); n=modes(k,2); typ=modes(k,3);
        ma=m*a; lam=ma^2+(n*pi)^2; Lam(k)=lam;
        sz=sin(n*pi*Zg); cz=cos(n*pi*Zg);
        if typ==0
            hx=cos(ma*Xg); dhx=-ma*sin(ma*Xg);
            E(k,:)=reshape(hx.*sz,1,[]);
            Ex(k,:)=reshape(dhx.*sz,1,[]);
            Ez(k,:)=reshape(n*pi*hx.*cz,1,[]);
            if m>0
                U(k,:)=reshape(-R*ma*n*pi/lam*sin(ma*Xg).*cz,1,[]);
                Wv(k,:)=reshape(R*ma^2/lam*cos(ma*Xg).*sz,1,[]);
            end
        else
            hx=sin(ma*Xg); dhx=ma*cos(ma*Xg);
            E(k,:)=reshape(hx.*sz,1,[]);
            Ex(k,:)=reshape(dhx.*sz,1,[]);
            Ez(k,:)=reshape(n*pi*hx.*cz,1,[]);
            U(k,:)=reshape(R*ma*n*pi/lam*cos(ma*Xg).*cz,1,[]);
            Wv(k,:)=reshape(R*ma^2/lam*sin(ma*Xg).*sz,1,[]);
        end
        if m==0, norms(k)=Lx/2; else, norms(k)=Lx/4; end
        if m==1 && n==1 && typ==0, idx11c=k; end
        if m==1 && n==1 && typ==1, idx11s=k; end
        id=3*(k-1)+(1:3);
        Lfull(id,id)=LmodalBase(Ra,a^2,H,J,M,m,n,A,gamma);
    end
    Proj=(E.*W2)./norms;
    Q=zeros(N,N,N);
    for i=1:N
        adv=(U(i,:).*Ex)+(Wv(i,:).*Ez);
        tmp=Proj*adv.';
        Q(:,i,:)=reshape(tmp,[N 1 N]);
    end
    Qflat=reshape(Q,N,N*N);
    model=struct('Mx',Mx,'Nz',Nz,'a',a,'H',H,'J',J,'M',M,'Ra',Ra,'A',A, ...
        'epsPor',epsPor,'gamma',gamma,'R',R,'modes',modes,'N',N,'Lx',Lx, ...
        'Lam',Lam,'Lfull',Lfull,'Q',Q,'Qflat',Qflat,'idx11c',idx11c,'idx11s',idx11s);
end

function Nv=nonlinearVector(model,y)
    X=reshape(y,3,model.N).'; th=X(:,1);
    % Vectorized quadratic contraction.  reshape(Q,N,N^2) uses the same
    % (i,j) ordering as kron(th,th); the product th(i)th(j) is symmetric.
    pairProducts=kron(th,th);
    Nth=-model.Qflat*pairProducts;
    Nmat=zeros(model.N,3); Nmat(:,1)=Nth;
    Nv=reshape(Nmat.',[],1);
end

function Nv=nonlinearVectorDirect(model,y)
    X=reshape(y,3,model.N).'; th=X(:,1);
    Nth=zeros(model.N,1);
    for i=1:model.N
        for j=1:model.N
            Nth=Nth-model.Q(:,i,j)*th(i)*th(j);
        end
    end
    Nmat=zeros(model.N,3); Nmat(:,1)=Nth;
    Nv=reshape(Nmat.',[],1);
end

function f=rhsModel(model,y)
    f=model.Lfull*y+nonlinearVector(model,y);
end

function Jm=modelJacobian(model,y)
    X=reshape(y,3,model.N).'; th=X(:,1);
    Jm=model.Lfull;
    D=zeros(model.N);
    for k=1:model.N
        for l=1:model.N
            s=0;
            for j=1:model.N
                s=s+(model.Q(k,l,j)+model.Q(k,j,l))*th(j);
            end
            D(k,l)=-s;
        end
    end
    for k=1:model.N
        for l=1:model.N
            Jm(3*(k-1)+1,3*(l-1)+1)=Jm(3*(k-1)+1,3*(l-1)+1)+D(k,l);
        end
    end
end

function [ts,Ys,Nus]=integrateSBDF2(model,X0,dt,T,saveDt)
    N3=3*model.N; nsteps=round(T/dt); saveEvery=max(1,round(saveDt/dt));
    I=eye(N3);
    D1=decomposition(I-dt*model.Lfull,'lu');
    D2=decomposition(3*I-2*dt*model.Lfull,'lu');
    y0=reshape(X0.',[],1); N0=nonlinearVector(model,y0);
    y1=D1\(y0+dt*N0); N1=nonlinearVector(model,y1);
    nsave=floor(nsteps/saveEvery)+1;
    ts=zeros(nsave,1); Ys=zeros(nsave,model.N,3); Nus=zeros(nsave,1);
    Ys(1,:,:)=reshape(X0,1,model.N,3); Nus(1)=compositeNu(model,X0); is=2;
    ym=y0; yn=y1; Nm=N0; Nn=N1;
    for k=1:nsteps-1
        yp=D2\(4*yn-ym+2*dt*(2*Nn-Nm));
        Np=nonlinearVector(model,yp);
        if mod(k+1,saveEvery)==0
            Xp=reshape(yp,3,model.N).';
            ts(is)=(k+1)*dt; Ys(is,:,:)=reshape(Xp,1,model.N,3);
            Nus(is)=compositeNu(model,Xp); is=is+1;
        end
        ym=yn; yn=yp; Nm=Nn; Nn=Np;
    end
    ts=ts(1:is-1); Ys=Ys(1:is-1,:,:); Nus=Nus(1:is-1);
end

function Nu=compositeNu(model,X)
    th=X(:,1); s=0;
    for k=1:model.N
        m=model.modes(k,1);
        if m>0, s=s+(m*model.a)^2/model.Lam(k)*th(k)^2; end
    end
    Nu=1+s/(4*(1+1/model.gamma));
end

function z=centerCoordinate(model,y,p)
    id=3*(model.idx11c-1)+(1:3);
    z=p'*y(id);
end

function [xg,Tg]=estimateOrbitGuess(ts,Ys,model,p,tWindow)
    keep=ts>=ts(end)-tWindow;
    t=ts(keep); Y=Ys(keep,:,:); z=zeros(numel(t),1);
    for k=1:numel(t)
        y=reshape(squeeze(Y(k,:,:)).',[],1);
        z(k)=centerCoordinate(model,y,p);
    end
    ph=unwrap(angle(z)); pp=polyfit(t,ph,1); om=abs(pp(1)); Tg=2*pi/om;
    last=t>=t(end)-Tg; ids=find(keep); ids=ids(last);
    zr=zeros(numel(ids),1);
    for k=1:numel(ids)
        y=reshape(squeeze(Ys(ids(k),:,:)).',[],1);
        zr(k)=real(centerCoordinate(model,y,p));
    end
    [~,jj]=max(zr);
    xg=reshape(squeeze(Ys(ids(jj),:,:)).',[],1);
end

function [x,T,Phi,res,it]=shootCorrect(model,x0,T0,maxit,tol)
    x=x0; T=T0; ref=x0; fref=rhsModel(model,ref); d=numel(x);
    for it=1:maxit
        [xT,Phi]=flowMonodromy(model,x,T);
        F=[xT-x; fref'*(x-ref)];
        res=norm(F);
        if res<tol, return; end
        J=[Phi-eye(d), rhsModel(model,xT); fref', 0];
        delta=J\(-F);
        x=x+delta(1:d); T=T+delta(end);
        if T<=0, T=abs(T); end
    end
    [xT,Phi]=flowMonodromy(model,x,T);
    res=norm([xT-x; fref'*(x-ref)]);
end

function [xT,Phi]=flowMonodromy(model,x,T)
    d=numel(x); Y0=[x; reshape(eye(d),[],1)];
    opts=odeset('RelTol',2e-9,'AbsTol',1e-11,'MaxStep',min(0.003,T/800));
    [~,Y]=ode45(@(t,Y) augmentedRhs(t,Y,model),[0 T],Y0,opts);
    yf=Y(end,:).'; xT=yf(1:d); Phi=reshape(yf(d+1:end),d,d);
end

function dY=augmentedRhs(~,Y,model)
    d=3*model.N; y=Y(1:d); Phi=reshape(Y(d+1:end),d,d);
    f=rhsModel(model,y); J=modelJacobian(model,y);
    dY=[f; reshape(J*Phi,[],1)];
end

function dg=orbitDiagnostics(model,x,T,p)
    tt=linspace(0,T,1001);
    opts=odeset('RelTol',2e-10,'AbsTol',1e-12,'MaxStep',0.003);
    [~,Y]=ode45(@(~,y) rhsModel(model,y),tt,x,opts);
    z=zeros(numel(tt),1); Nu=zeros(numel(tt),1);
    th=zeros(numel(tt),1); sf=zeros(numel(tt),1); d2=zeros(numel(tt),1);
    weights=0.25*ones(model.N,1); weights(model.modes(:,1)==0)=0.5;
    for k=1:numel(tt)
        z(k)=centerCoordinate(model,Y(k,:).',p);
        X=reshape(Y(k,:),3,model.N).';
        Nu(k)=compositeNu(model,X);
        th(k)=X(model.idx11c,1); sf(k)=X(model.idx11c,2);
        d=X(:,1)-X(:,2); d2(k)=sum(weights.*d.^2);
    end
    omega=2*pi/T;
    th0=th-trapz(tt,th)/T; sf0=sf-trapz(tt,sf)/T;
    e=exp(-1i*omega*tt(:));
    cth=trapz(tt,th0.*e)/T; csf=trapz(tt,sf0.*e)/T;
    dp=angle(csf)-angle(cth); dp=mod(dp+pi,2*pi)-pi;
    dg=struct('r',trapz(tt,abs(z))/T,'Nu',trapz(tt,Nu)/T, ...
        'phaseGap',abs(dp),'DT',sqrt(trapz(tt,d2)/T), ...
        'NuHalfAmp',0.5*(max(Nu)-min(Nu)));
end

function dg=diagnosticsPeriodic(ts,Ys,Nus,model,p,tmin)
    id=find(ts>=tmin); t=ts(id); Y=Ys(id,:,:); Nu=Nus(id);
    zc=zeros(numel(id),1);
    for k=1:numel(id)
        xcrit=squeeze(Y(k,model.idx11c,:)); zc(k)=p'*xcrit;
    end
    ph=unwrap(angle(zc)); pp=polyfit(t,ph,1); omega=pp(1);
    r=trapz(t,abs(zc))/(t(end)-t(1));
    th=squeeze(Y(:,model.idx11c,1)); sf=squeeze(Y(:,model.idx11c,2));
    th0=th-trapz(t,th)/(t(end)-t(1)); sf0=sf-trapz(t,sf)/(t(end)-t(1));
    e=exp(-1i*omega*t);
    cth=trapz(t,th0.*e)/(t(end)-t(1)); csf=trapz(t,sf0.*e)/(t(end)-t(1));
    dp=angle(csf)-angle(cth); dp=mod(dp+pi,2*pi)-pi;
    weights=0.25*ones(model.N,1); weights(model.modes(:,1)==0)=0.5;
    d2=zeros(numel(id),1);
    for k=1:numel(id)
        X=squeeze(Y(k,:,:)); d=X(:,1)-X(:,2); d2(k)=sum(weights.*d.^2);
    end
    dg=struct('r',r,'omega',omega,'Nu',trapz(t,Nu)/(t(end)-t(1)), ...
        'phaseGap',abs(dp),'DT',sqrt(trapz(t,d2)/(t(end)-t(1))), ...
        'NuHalfAmp',0.5*(max(Nu)-min(Nu)));
end

function yFull=embedCosState(modelC,modelO,xC)
    XC=reshape(xC,3,modelC.N).'; XO=zeros(modelO.N,3);
    for kc=1:modelC.N
        mode=modelC.modes(kc,:);
        kf=find(modelO.modes(:,1)==mode(1) & modelO.modes(:,2)==mode(2) & modelO.modes(:,3)==0,1);
        XO(kf,:)=XC(kc,:);
    end
    yFull=reshape(XO.',[],1);
end

function [muMax,exponent,muNeutral]=transverseFloquet(modelO,y0,T)
    sineModes=find(modelO.modes(:,3)==1);
    idx=[];
    for k=1:numel(sineModes), idx=[idx,3*(sineModes(k)-1)+(1:3)]; end %#ok<AGROW>
    ds=numel(idx); Y0=[y0; reshape(eye(ds),[],1)];
    opts=odeset('RelTol',2e-9,'AbsTol',1e-11,'MaxStep',min(0.003,T/800));
    [~,Y]=ode45(@(t,Y) transverseAugRhs(t,Y,modelO,idx),[0 T],Y0,opts);
    yf=Y(end,:).'; Phi=reshape(yf(numel(y0)+1:end),ds,ds);
    mu=eig(Phi);
    [~,i0]=min(abs(mu-1)); muNeutral=abs(mu(i0)); mu(i0)=[];
    muMax=max(abs(mu)); exponent=log(muMax)/T;
end

function dY=transverseAugRhs(~,Y,model,idx)
    d=3*model.N; ds=numel(idx);
    y=Y(1:d); Phi=reshape(Y(d+1:end),ds,ds);
    f=rhsModel(model,y); J=modelJacobian(model,y); Js=J(idx,idx);
    dY=[f; reshape(Js*Phi,[],1)];
end

function dg=travelingDiagnostics(ts,Ys,Nus,model,tmin)
    keep=ts>=tmin; t=ts(keep); Y=Ys(keep,:,:); Nu=Nus(keep);
    anti=zeros(numel(t),1); phase=zeros(numel(t),1);
    sineModes=find(model.modes(:,3)==1);
    for k=1:numel(t)
        X=squeeze(Y(k,:,:));
        anti(k)=norm(X(sineModes,:),'fro');
        phase(k)=atan2(X(model.idx11s,1),X(model.idx11c,1));
    end
    phase=unwrap(phase); pp=polyfit(t,phase,1);
    dg.speed=pp(1)/model.a;
    dg.Nu=trapz(t,Nu)/(t(end)-t(1));
    dg.anti=fullAntisymmetricSeries(Ys,model);
    % Growth estimate from the complete time series when amplitudes are small.
    aa=dg.anti; mask=aa>2e-4 & aa<2e-2;
    if sum(mask)>10
        pp2=polyfit(ts(mask),log(aa(mask)),1); dg.growth=pp2(1);
    else
        dg.growth=NaN;
    end
end

function anti=fullAntisymmetricSeries(Ys,model)
    sineModes=find(model.modes(:,3)==1);
    anti=zeros(size(Ys,1),1);
    for k=1:size(Ys,1)
        X=squeeze(Ys(k,:,:));
        anti(k)=norm(X(sineModes,:),'fro');
    end
end

function [x,z,theta,phi]=reconstructFields(model,Xcoef,Nx,Nz)
    x=linspace(0,model.Lx,Nx); z=linspace(0,1,Nz);
    [XX,ZZ]=meshgrid(x,z); theta=zeros(size(XX)); phi=zeros(size(XX));
    for k=1:model.N
        m=model.modes(k,1); n=model.modes(k,2); typ=model.modes(k,3);
        if typ==0, hx=cos(m*model.a*XX); else, hx=sin(m*model.a*XX); end
        b=hx.*sin(n*pi*ZZ);
        theta=theta+Xcoef(k,1)*b; phi=phi+Xcoef(k,2)*b;
    end
end

function inter=steadyHopfInteraction(H,J,A,gamma)
    [RaStar,a2s]=critStat(H,gamma);
    ks=sqrt(a2s); ms=47; mh=46; a0=ks/ms; base2=a0^2;
    a2h=(mh*a0)^2;
    Mstar=fzero(@(M) hopfRaAtA2(a2h,H,J,M,A,gamma)-RaStar,[1.245 1.25], ...
        optimset('TolX',1e-12,'Display','off'));
    Ls=LmodalBase(RaStar,base2,H,J,Mstar,ms,1,A,gamma);
    Lh=LmodalBase(RaStar,base2,H,J,Mstar,mh,1,A,gamma);
    [Vs,Ds]=eig(complex(Ls)); es=diag(Ds); [~,is]=min(abs(es));
    s=Vs(:,is); s=s/norm(s); s=s*exp(-1i*angle(s(1))); s=real(s);
    [Ps,Dps]=eig(Ls'); epsv=diag(Dps); [~,ips]=min(abs(epsv));
    ps=Ps(:,ips); ps=real(ps/(ps'*s));
    [Vh,Dh]=eig(complex(Lh)); eh=diag(Dh); pos=find(imag(eh)>0);
    [~,jj]=min(abs(real(eh(pos)))); ih=pos(jj); q=Vh(:,ih); q=q/norm(q); q=q*exp(-1i*angle(q(1)));
    omega=imag(eh(ih));
    [Ph,Dph]=eig(Lh'); eph=diag(Dph); [~,iph]=min(abs(eph+1i*omega));
    p=Ph(:,iph); p=p/conj(p'*q);

    NX=256; NZ=28;
    bss=BcoefCos(ms,1,ms,1,0,2,a0,RaStar,NX,NZ)*s(1)^2;
    bqq=BcoefCos(mh,1,mh,1,0,2,a0,RaStar,NX,NZ)*q(1)^2;
    bqb=BcoefCos(mh,1,mh,1,0,2,a0,RaStar,NX,NZ)*q(1)*conj(q(1));
    L02=LmodalBase(RaStar,base2,H,J,Mstar,0,2,A,gamma);
    hss=-L02\[bss;0;0];
    hqq=(2i*omega*eye(3)-L02)\[bqq;0;0];
    hqb=-L02\[bqb;0;0];

    mts=[1 93]; hsq=cell(2,1);
    for r=1:2
        mt=mts(r);
        bsq=BcoefCos(ms,1,mh,1,mt,2,a0,RaStar,NX,NZ)*s(1)*q(1);
        Lx=LmodalBase(RaStar,base2,H,J,Mstar,mt,2,A,gamma);
        hsq{r}=(1i*omega*eye(3)-Lx)\[bsq;0;0];
    end

    bS0=BcoefCos(ms,1,0,2,ms,1,a0,RaStar,NX,NZ);
    a_s=0.5*real(ps(1)*bS0*s(1)*hss(1));

    bH0=BcoefCos(mh,1,0,2,mh,1,a0,RaStar,NX,NZ);
    c_h=p'*[bH0*(q(1)*hqb(1)+0.5*conj(q(1))*hqq(1));0;0];

    fb_hs=0.5*bH0*q(1)*hss(1);
    for r=1:2
        fb_hs=fb_hs+BcoefCos(ms,1,mts(r),2,mh,1,a0,RaStar,NX,NZ)*s(1)*hsq{r}(1);
    end
    c_hs=p'*[fb_hs;0;0];

    fb_sh=bS0*s(1)*hqb(1);
    for r=1:2
        bc=BcoefCos(mh,1,mts(r),2,ms,1,a0,RaStar,NX,NZ);
        fb_sh=fb_sh+bc*(q(1)*conj(hsq{r}(1))+conj(q(1))*hsq{r}(1));
    end
    a_sh=real(ps(1)*fb_sh);

    LraS=zeros(3); LamS=(ms*a0)^2+pi^2; LraS(1,1)=(ms*a0)^2/LamS;
    kappa_s=real(ps'*(LraS*s));
    LraH=zeros(3); LamH=(mh*a0)^2+pi^2; LraH(1,1)=(mh*a0)^2/LamH;
    kh=p'*(LraH*q);
    dM=1e-6;
    dLM=(LmodalBase(RaStar,base2,H,J,Mstar+dM,mh,1,A,gamma)- ...
         LmodalBase(RaStar,base2,H,J,Mstar-dM,mh,1,A,gamma))/(2*dM);
    eta=p'*(dLM*q);

    detMix=a_s*real(c_h)-a_sh*real(c_hs);
    hopfExist=-real(kh)/real(eta);
    lower=(real(c_h)*kappa_s/a_sh-real(kh))/real(eta);
    upper=-(real(kh)-real(c_hs)*kappa_s/a_s)/real(eta);

    inter=struct('Ra_star',RaStar,'M_star',Mstar,'omega_star',omega,'a0',a0, ...
        'Lx',2*pi/a0,'m_stationary',ms,'m_hopf',mh,'k_stationary',ks, ...
        'k_hopf',mh*a0,'a_s',a_s,'a_sh',a_sh,'c_h_real',real(c_h), ...
        'c_h_imag',imag(c_h),'c_hs_real',real(c_hs),'c_hs_imag',imag(c_hs), ...
        'determinant',detMix,'kappa_s',kappa_s,'kappa_h_real',real(kh), ...
        'kappa_h_imag',imag(kh),'eta_h_real',real(eta),'eta_h_imag',imag(eta), ...
        'hopf_exist_slope',hopfExist,'lower_slope',lower,'upper_slope',upper);
end

function c=BcoefCos(mi,ni,mj,nj,mt,nt,a0,Ra,Nx,Nz)
    Lx=2*pi/a0; x=(0:Nx-1)'*Lx/Nx; [z,wz]=gaussLegendre(Nz,0,1);
    [X,Z]=ndgrid(x,z); W=(Lx/Nx)*ones(Nx,1)*wz.';
    [ei,exi,ezi,ui,wi]=basisCos(mi,ni,a0,Ra,X,Z);
    [ej,exj,ezj,uj,wj]=basisCos(mj,nj,a0,Ra,X,Z);
    [et,~,~,~,~]=basisCos(mt,nt,a0,Ra,X,Z);
    b=-(ui.*exj+wi.*ezj+uj.*exi+wj.*ezi);
    if mt==0, nrm=Lx/2; else, nrm=Lx/4; end
    c=sum(W.*b.*et,'all')/nrm;
end

function [e,ex,ez,u,w]=basisCos(m,n,a0,Ra,X,Z)
    k=m*a0; Lam=k^2+(n*pi)^2;
    e=cos(k*X).*sin(n*pi*Z);
    ex=-k*sin(k*X).*sin(n*pi*Z);
    ez=n*pi*cos(k*X).*cos(n*pi*Z);
    if m==0
        u=zeros(size(X)); w=zeros(size(X));
    else
        u=-sqrt(Ra)*k*n*pi/Lam*sin(k*X).*cos(n*pi*Z);
        w=sqrt(Ra)*k^2/Lam*cos(k*X).*sin(n*pi*Z);
    end
end

function [RaE,a2]=energyThreshold(H,J,M,gamma)
    f=@(x) ((pi^2+x).^2./x)+H*(pi^2+x).^2./ ...
        (x.*(pi^2+x+gamma*H*(M/J)));
    [a2,RaE]=fminbnd(f,2,30,optimset('TolX',1e-11,'Display','off'));
end

function [x,w]=gaussLegendre(n,a,b)
    k=(1:n-1)'; beta=k./sqrt(4*k.^2-1);
    T=diag(beta,1)+diag(beta,-1);
    [V,D]=eig(T); [x,ix]=sort(diag(D)); V=V(:,ix);
    w=2*(V(1,:).^2).';
    x=(a+b)/2+(b-a)/2*x; w=(b-a)/2*w;
end

function C=publicationColors()
    C=[0.85 0.10 0.10; ...
       0.00 0.32 0.82; ...
       0.00 0.00 0.00; ...
       0.55 0.20 0.72; ...
       0.00 0.50 0.30];
end

function styleAxes(ax)
    set(ax,'FontSize',14,'TickLabelInterpreter','latex','LineWidth',1.0, ...
        'Box','on','Layer','top');
    ax.XLabel.Interpreter='latex'; ax.YLabel.Interpreter='latex';
    ax.XLabel.FontSize=14; ax.YLabel.FontSize=14;
end

function lgd=placeLegend(ax,loc)
    if nargin<2 || isempty(loc), loc='northeast'; end
    loc=strrep(loc,'outside','');
    lgd=legend(ax,'Location',loc,'Interpreter','latex','FontSize',12, ...
        'Box','off','Color','none','EdgeColor','none');
    try, lgd.ItemTokenSize=[24 12]; catch, end
end


function cmap=academicSequentialMap(n)
    % Conservative Q1-journal sequential palette for one-signed scalar fields.
    % Dark navy denotes the lower end of the range and pale blue the upper end.
    % The monotone lightness progression remains interpretable in grayscale.
    if nargin<1, n=256; end
    anchors=[0.0314 0.1882 0.4196; ... % deep navy
             0.1294 0.4431 0.7098; ... % academic blue
             0.2588 0.5725 0.7765; ... % medium blue
             0.4196 0.6824 0.8392; ... % light blue
             0.6196 0.7922 0.8824; ... % pale blue
             0.8706 0.9216 0.9686];    % very pale blue
    xa=linspace(0,1,size(anchors,1)); x=linspace(0,1,n); cmap=zeros(n,3);
    for j=1:3, cmap(:,j)=interp1(xa,anchors(:,j),x,'pchip'); end
    cmap=max(0,min(1,cmap));
end

function cmap=academicDivergingMap(n)
    % Conservative blue--neutral--red diverging palette for signed fields.
    % Negative values are blue, zero is near white, positive values are red.
    % Symmetric colour limits must be used with this map.
    if nargin<1, n=256; end
    anchors=[0.0196 0.1882 0.3804; ... % deep blue
             0.1294 0.4000 0.6745; ... % blue
             0.5725 0.7725 0.8706; ... % pale blue
             0.9686 0.9686 0.9686; ... % neutral near white
             0.9569 0.6471 0.5098; ... % pale red
             0.8392 0.3765 0.3020; ... % red
             0.4039 0.0000 0.1216];    % deep red
    xa=linspace(0,1,size(anchors,1)); x=linspace(0,1,n); cmap=zeros(n,3);
    for j=1:3, cmap(:,j)=interp1(xa,anchors(:,j),x,'pchip'); end
    cmap=max(0,min(1,cmap));
end


function writeFinalLatexTables(root,Tval,Tmap,Tcont,Ttrans,Ttw,Tsens,Tsp,Ttime)
    % Core manuscript tables. The independent-grid/refined-3D table is
    % produced by export_independent_validation_HFF.m and inserted at the
    % marker below so tables_generated.tex matches the final manuscript.
    fid=fopen(fullfile(root,'tables_core_generated.tex'),'w');
    if fid<0, error('Unable to create tables_core_generated.tex'); end
    cleaner=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%% Auto-generated core tables by run_all_results_HFF_final.m.\n');
    fprintf(fid,'%% The independent-grid/refined-3D table is inserted by export_independent_validation_HFF.m.\n\n');

    %% Table I: parent-theory validation
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Independent recovery of selected stationary--oscillatory transition data from \\citet{Straughan2026}. \\textit{Source: Author''s calculations.}}\n');
    fprintf(fid,'\\label{tab:validation}\n\\small\n\\begin{tabular}{rrrrrrr}\n\\toprule\n');
    fprintf(fid,'$H$ & $J$ & $(M/J)_t$ ref. & $(M/J)_t$ calc. & difference (\\%%) & $\\mathrm{Ra}_{\\mathrm{stat}}$ ref. & $\\mathrm{Ra}_{\\mathrm{stat}}$ calc.\\\\\n\\midrule\n');
    for k=1:height(Tval)
        fprintf(fid,'%g & %.2f & %.4f & %.4f & %.3f & %.4f & %.4f\\\\\n', ...
            Tval.H(k),Tval.J(k),Tval.parent_ratio(k),Tval.computed_ratio(k), ...
            Tval.ratio_error_pct(k),Tval.parent_Rastat(k),Tval.computed_Rastat(k));
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');

    %% Former Table II removed: Hopf-survey ranges are reported compactly in the Introduction.

    %% Table II: periodic continuation
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Selected periodic-orbit continuation data for the benchmark case. \\textit{Source: Author''s calculations.}}\n');
    fprintf(fid,'\\label{tab:continuation}\n\\small\n\\begin{tabular}{rrrrrr}\n\\toprule\n');
    fprintf(fid,'$\\mathrm{Ra}-\\mathrm{Ra}_c$ & $r$ & $\\omega$ & $\\overline{Nu}$ & $|\\mu_F|_{\\max}$ & shooting residual\\\\\n\\midrule\n');
    for k=1:height(Tcont)
        fprintf(fid,'%.2f & %.6f & %.6f & %.6f & %.6f & $%.1e$\\\\\n', ...
            Tcont.deltaRa(k),Tcont.r_shoot(k),Tcont.omega_shoot(k),Tcont.Nu_mean(k), ...
            Tcont.max_nontrivial_Floquet(k),Tcont.shooting_residual(k));
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');

    %% Table III: transverse Floquet and traveling-wave diagnostics
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Dominant transverse Floquet multiplier and traveling-wave diagnostics in the extended two-dimensional calculation. \\textit{Source: Author''s calculations.}}\n');
    fprintf(fid,'\\label{tab:transverse}\n\\small\n\\begin{tabular}{rrrrr}\n\\toprule\n');
    fprintf(fid,'$\\mathrm{Ra}-\\mathrm{Ra}_c$ & $|\\mu_\\perp|$ & transverse exponent & $|c_{\\mathrm{TW}}|$ & $\\overline{Nu}_{\\mathrm{TW}}$\\\\\n\\midrule\n');
    for k=1:height(Ttrans)
        fprintf(fid,'%.2f & %.6f & %.6f & %.6f & %.6f\\\\\n', ...
            Ttrans.deltaRa(k),Ttrans.max_transverse_multiplier_abs(k),Ttrans.transverse_exponent(k), ...
            Ttw.wave_speed_abs(k),Ttw.Nu_traveling(k));
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');

    %% Table IV: targeted nonlinear parameter sensitivity
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Targeted nonlinear sensitivity at $\\mathrm{Ra}-\\mathrm{Ra}_c=1$. The benchmark and two additional constitutive states are evaluated with the same two-dimensional Galerkin formulation. Positive $\\sigma_\\perp$ denotes transverse growth of the standing periodic orbit. \\textit{Source: Author''s calculations.}}\n');
    fprintf(fid,'\\label{tab:sensitivity}\n\\small\n\\begin{tabular}{rrrrrrrr}\n\\toprule\n');
    fprintf(fid,'$H$ & $J$ & $M$ & $M/M_t$ & $\\mathrm{Ra}_c$ & $\\sigma_\\perp$ & $|c_{\\rm TW}|$ & $\\overline{Nu}_{\\rm SW}\\;/\\;\\overline{Nu}_{\\rm TW}$\\\\\n\\midrule\n');
    for k=1:height(Tsens)
        fprintf(fid,'%g & %.1f & %.1f & %.3f & %.5f & %.5f & %.6f & %.6f / %.6f\\\\\n', ...
            Tsens.H(k),Tsens.J(k),Tsens.M(k),Tsens.M_over_Mt(k),Tsens.Rac(k), ...
            Tsens.transverse_exponent(k),Tsens.wave_speed_abs(k),Tsens.Nu_standing(k),Tsens.Nu_traveling(k));
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');

    fprintf(fid,'%% __GRID_TABLE_INSERTION_POINT__\n\n');

    %% Table VI after independent-grid table insertion: convergence
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Spatial and temporal convergence at $\\mathrm{Ra}-\\mathrm{Ra}_c=0.2$. Spatial refinement uses corrected periodic shooting, including a demanding $(6,8)$ check; temporal refinement uses $(M_x,N_z)=(5,6)$ and successively halved time steps. \\textit{Source: Author''s calculations.}}\n');
    fprintf(fid,'\\label{tab:convergence}\n\\small\n\\begin{tabular}{llrrrr}\n\\toprule\n');
    fprintf(fid,'study & resolution & period or $\\Delta t$ & $r$ & $\\omega$ & $\\overline{Nu}$\\\\\n\\midrule\n');
    for k=1:height(Tsp)
        fprintf(fid,'spatial & $(%d,%d)$ & %.9f & %.9f & %.9f & %.9f\\\\\n', ...
            Tsp.Mx(k),Tsp.Nz(k),Tsp.period(k),Tsp.r(k),Tsp.omega(k),Tsp.Nu(k));
    end
    fprintf(fid,'\\midrule\n');
    for k=1:height(Ttime)
        fprintf(fid,'time & $(5,6)$ & %.4f & %.8f & %.8f & %.8f\\\\\n', ...
            Ttime.dt(k),Ttime.r(k),Ttime.omega(k),Ttime.Nu(k));
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end

function saveFive(fig,root,name,rasterHeavy)
    if nargin<4, rasterHeavy=false; end
    drawnow;
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
    close(fig);
end
