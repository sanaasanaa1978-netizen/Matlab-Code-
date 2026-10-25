%% publication_postprocess_HFF.m
% Publication-strength post-processing for the HFF reproducibility package.
% Run after run_all_results_HFF_final.m and independent_grid_validation_HFF.m.
% Outputs:
%   independent_validation/richardson_convergence_analysis.csv
%   numerical_uncertainty_budget.csv
%   transition_ratio_resolution_matlab.csv
%
% Numerical uncertainties below are deterministic discretization/solver
% estimates, not statistical confidence intervals.

clear; clc;
root=fileparts(mfilename('fullpath'));
outdir=fullfile(root,'independent_validation');
G=readtable(fullfile(outdir,'independent_grid_convergence_matlab.csv'));
GR=readtable(fullfile(outdir,'independent_grid_convergence_matlab_with_reference.csv'));

% The grids satisfy Nx=2*(Nz+1), so h=1/(Nz+1) scales both dx and dz.
h=1./(G.Nz+1);
qNames={'wave_speed_abs','omega_traveling','Nu_traveling'};
obs=zeros(3,1); finest=zeros(3,1); extrap=zeros(3,1); eabs=zeros(3,1); epct=zeros(3,1);
gref=zeros(3,1); gabs=zeros(3,1); gpct=zeros(3,1);
for iq=1:3
    f=G.(qNames{iq});
    obs(iq)=observedOrderUnequal(h(1:3),f(1:3));
    finest(iq)=f(4);
    extrap(iq)=f(4)+(f(4)-f(3))/((h(3)/h(4))^2-1); % scheme-order p=2
    eabs(iq)=abs(extrap(iq)-finest(iq));
    epct(iq)=100*eabs(iq)/abs(extrap(iq));
    gref(iq)=GR.(qNames{iq})(end);
    gabs(iq)=abs(finest(iq)-gref(iq));
    gpct(iq)=100*gabs(iq)/abs(gref(iq));
end
TR=table(string(qNames(:)),obs,finest,extrap,eabs,epct,gref,gabs,gpct, ...
    'VariableNames',{'quantity','observed_order_32_48_64','finest_value_80x39', ...
    'second_order_extrapolated_value_64_80','finest_to_extrap_abs', ...
    'finest_to_extrap_pct','galerkin_reference','finest_to_galerkin_abs','finest_to_galerkin_pct'});
writetable(TR,fullfile(outdir,'richardson_convergence_analysis.csv'));
writetable(TR,fullfile(root,'richardson_convergence_analysis.csv'));

% Direct check of the one parent-study transition value with the largest mismatch.
H=10; J=0.5; A=0.13164; gamma=0.3/0.7;
[Rstat,aStat2]=critStatPP(H,gamma); %#ok<ASGLU>
rPub=2.501; MPub=J*rPub;
[RaPub,a2Pub,omPub]=critHopfPP(H,J,MPub,A,gamma);
rStar=fzero(@(r) transitionResidualPP(r,H,J,A,gamma,Rstat),[2.49 2.50],optimset('TolX',1e-12,'Display','off'));
MStar=J*rStar; [RaStar,a2Star,omStar]=critHopfPP(H,J,MStar,A,gamma);
leadPub=leadingRealPP(Rstat,a2Pub,H,J,MPub,A,gamma);
leadStar=leadingRealPP(Rstat,a2Star,H,J,MStar,A,gamma);
caseName=["published_ratio";"direct_recomputed_crossing"];
TH=H*ones(2,1); TJ=J*ones(2,1); ratio=[rPub;rStar]; MM=[MPub;MStar];
RaS=Rstat*ones(2,1); RaO=[RaPub;RaStar]; dRa=RaO-RaS; a2=[a2Pub;a2Star];
omega=[omPub;omStar]; lead=[leadPub;leadStar];
TT=table(caseName,TH,TJ,ratio,MM,RaS,RaO,dRa,a2,omega,lead,'VariableNames', ...
 {'case','H','J','M_over_J','M','Ra_stat','Ra_osc_min','Ra_osc_minus_Ra_stat','a2_osc','omega_osc','leading_real_part_at_Ra_stat'});
writetable(TT,fullfile(root,'transition_ratio_resolution_matlab.csv'));

% Compact deterministic numerical uncertainty budget.
V=readtable(fullfile(root,'matlab_run_verification.csv'));
TC=readtable(fullfile(root,'transverse_convergence_matlab.csv'));
S=readtable(fullfile(root,'limited_3D_spanwise_stability_matlab.csv'));
idx190=find(abs(S.ky_over_a-0.190)<1e-12,1); idx192=find(abs(S.ky_over_a-0.192)<1e-12,1);
s1=S.growth_exponent_48x23(idx190); s2=S.growth_exponent_48x23(idx192);
kcross=0.190-s1*(0.002)/(s2-s1);
mu=TC.max_transverse_multiplier_abs(end);
muUnc=abs(mu-TC.max_transverse_multiplier_abs(end-1));
getv=@(nm) V.value(strcmp(V.check,nm));
gett=@(nm) V.tolerance(strcmp(V.check,nm));
rowQ=["critical Rayleigh number Ra_c";"Hopf frequency omega_c";"first Lyapunov coefficient ell1"; ...
      "traveling-wave speed |c| at deltaRa=1 (independent grid)"; ...
      "traveling Nusselt number at deltaRa=1 (independent grid)"; ...
      "dominant transverse Floquet multiplier |mu_perp| at deltaRa=0.2"; ...
      "spanwise neutral crossover k_y/a"];
rv=[getv('benchmark_Rac');getv('benchmark_omega');getv('benchmark_ell1');G.wave_speed_abs(end);G.Nu_traveling(end);mu;kcross];
ua=[gett('benchmark_Rac');gett('benchmark_omega');gett('benchmark_ell1'); ...
    max(TR.finest_to_extrap_abs(1),TR.finest_to_galerkin_abs(1)); ...
    max(TR.finest_to_extrap_abs(3),TR.finest_to_galerkin_abs(3));muUnc;0.001];
ur=100*ua./abs(rv);
basis=["conservative benchmark/reoptimization verification envelope; 3x3 linear eigenproblem"; ...
       "conservative benchmark/reoptimization verification envelope; 3x3 linear eigenproblem"; ...
       "conservative normal-form algebra verification envelope"; ...
       "max of finest-grid to p=2 Richardson extrapolation and Galerkin cross-discretization difference"; ...
       "max of finest-grid to p=2 Richardson extrapolation and Galerkin cross-discretization difference"; ...
       "absolute change from (4,5) to (6,8) Galerkin resolution"; ...
       "half-width of resolved sign-change bracket [0.190,0.192]"];
TU=table(rowQ,rv,ua,ur,basis,'VariableNames',{'quantity','reported_value','estimated_abs_uncertainty','estimated_rel_uncertainty_pct','basis'});
writetable(TU,fullfile(root,'numerical_uncertainty_budget_matlab.csv'));

fprintf('Publication post-processing complete.\n');
fprintf('Observed order (32/48/64): speed %.4f, Nu %.4f.\n',TR.observed_order_32_48_64(1),TR.observed_order_32_48_64(3));
fprintf('Direct transition: published ratio %.6f -> dRa=%+.6g; recomputed ratio %.9f -> dRa=%+.3e.\n', ...
    rPub,RaPub-Rstat,rStar,RaStar-Rstat);

function p=observedOrderUnequal(h,f)
    R=(f(1)-f(2))/(f(2)-f(3));
    fun=@(x) (h(1)^x-h(2)^x)/(h(2)^x-h(3)^x)-R;
    p=fzero(fun,[0.05 10]);
end
function [Rstat,a2]=critStatPP(H,gamma)
    f=@(x) ((pi^2+x).^2./x)+H*(pi^2+x).^2./(x.*(pi^2+x+gamma*H));
    [a2,Rstat]=fminbnd(f,2,30,optimset('TolX',1e-11,'Display','off'));
end
function Ra=hopfRaAtA2PP(a2,H,J,M,A,gamma)
    Lam=a2+pi^2; alpha=(1/M)*(1-J/M);
    d=-((J/M)*Lam+gamma*H)/A; e=-alpha*Lam/A; ff=-1/M; b=gamma*H/A;
    u=d+ff; w=d*ff-e; v=w-H*b;
    rr=roots([-u,H*b-u^2,-(u*v+H*b*ff)]); cand=[];
    for k=1:numel(rr)
        if abs(imag(rr(k)))<1e-9
            x=real(rr(k)); Rt=(x+Lam+H)*Lam/a2; c2=-(x+d+ff); c1=x*(d+ff)+d*ff-e-H*b;
            if Rt>0 && c2>0 && c1>0, cand(end+1)=Rt; end %#ok<AGROW>
        end
    end
    if isempty(cand), Ra=NaN; else, Ra=min(cand); end
end
function [Ra,a2,omega]=critHopfPP(H,J,M,A,gamma)
    ag=linspace(2,30,70); rg=nan(size(ag));
    for k=1:numel(ag), rg(k)=hopfRaAtA2PP(ag(k),H,J,M,A,gamma); end
    ok=find(isfinite(rg)); [~,jj]=min(rg(ok)); x0=ag(ok(jj)); lo=max(0.2,x0-2); hi=x0+2;
    obj=@(x) finiteObjPP(x,H,J,M,A,gamma);
    [a2,Ra]=fminbnd(obj,lo,hi,optimset('TolX',1e-10,'Display','off'));
    L=LPP(Ra,a2,H,J,M,A,gamma); c=poly(L); omega=sqrt(c(3));
end
function y=finiteObjPP(x,H,J,M,A,gamma)
    y=hopfRaAtA2PP(x,H,J,M,A,gamma); if ~isfinite(y), y=1e12; end
end
function g=transitionResidualPP(r,H,J,A,gamma,Rstat)
    [Ro,~,~]=critHopfPP(H,J,J*r,A,gamma); g=Ro-Rstat;
end
function L=LPP(Ra,a2,H,J,M,A,gamma)
    Lam=a2+pi^2; alpha=(1/M)*(1-J/M); buoy=Ra*a2/Lam;
    L=[buoy-Lam-H,H,0;gamma*H/A,-((J/M)*Lam+gamma*H)/A,-alpha*Lam/A;0,1,-1/M];
end
function g=leadingRealPP(Ra,a2,H,J,M,A,gamma)
    g=max(real(eig(LPP(Ra,a2,H,J,M,A,gamma))));
end
