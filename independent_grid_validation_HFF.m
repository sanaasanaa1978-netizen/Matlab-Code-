%% independent_grid_validation_HFF.m
% Independent grid-based PDE validation of the standing-to-traveling-wave
% transition. This implementation does NOT reuse the Galerkin quadratic
% interaction tensor. It uses:
%   - periodic FFT differentiation in x,
%   - centered second-order finite differences in z,
%   - diagonalization of the discrete Dirichlet Laplacian by a sine matrix,
%   - a discrete Poisson solve for the Darcy velocity,
%   - pointwise nonlinear advection,
%   - IMEX SBDF2 time integration.
%
% This MATLAB implementation is independent of the Galerkin interaction
% tensor used by the primary workflow. If an external cross-implementation
% reference named independent_grid_convergence.csv is placed in the output
% directory, the script also reports MATLAB-versus-reference differences;
% that optional file is not required for the manuscript reproduction.

clear; clc;
root=fileparts(mfilename('fullpath'));
outdir=fullfile(root,'independent_validation');
if ~exist(outdir,'dir'), mkdir(outdir); end

% Base four-grid sequence used in the manuscript. Set runExtraRefinement=true
% to append a 96 x 47 confirmation without changing the
% published four-grid data or the Richardson analysis.
runExtraRefinement=false;
cfg=[32 15 0.0100 160; ...
     48 23 0.0075 160; ...
     64 31 0.0050 180; ...
     80 39 0.0050 120];
if runExtraRefinement
    cfg=[cfg; 96 47 0.00375 120]; %#ok<AGROW>
end
rows=zeros(size(cfg,1),10);
allHigh=[]; spanBase32=[]; spanBase48=[]; spanBase64=[]; spanBase80=[];
for ic=1:size(cfg,1)
    Nx=cfg(ic,1); Nz=cfg(ic,2); dt=cfg(ic,3); Ttravel=cfg(ic,4);
    fprintf('Independent grid case %dx%d, dt=%g ...\n',Nx,Nz,dt);
    S=makeSolver(Nx,Nz,dt,1.0);
    F=initialField(S,0.70,0);
    [F,ts,rec]=advanceGrid(S,F,80,true,0.04,false);

    [X,Z]=ndgrid(S.x,S.z);
    for j=1:3
        F(:,:,j)=F(:,:,j)+1e-4*real(S.q(j))*sin(S.a*X).*sin(pi*Z);
    end
    saveSnapshots=(Nx==64 && Nz==31);
    [F2,tt,rr,snaps,snapTimes]=advanceGrid(S,F,Ttravel,false,0.04,saveSnapshots);

    % Standing frequency from a least-squares sinusoidal fit after transient.
    % This avoids the harmonic/secondary-maximum sensitivity of the earlier
    % local-maximum picker used in an earlier verification run.
    sig=rec(:,3).*cos(rec(:,2));  % amplitude * cos(FFT phase)
    use=ts>30; tm=ts(use); sm=sig(use)-mean(sig(use));
    omegaSW=fitOscillationFrequency(tm,sm,1.0,1.6);

    late=tt>(Ttravel-60);
    ph=-unwrap(rr(:,2)); pp=polyfit(tt(late),ph(late),1);
    cSigned=pp(1)/S.a; cTW=abs(cSigned); omegaTW=abs(pp(1));
    NuSW=mean(rec(use,4)); NuTW=mean(rr(late,4));
    [ct,st]=critProjection(S,F2,1); thetaAmp=hypot(ct,st);
    [cp,sp]=critProjection(S,F2,2);
    dphi=abs(wrapPi(atan2(sp,cp)-atan2(st,ct)));

    rows(ic,:)=[Nx,Nz,dt,omegaSW,NuSW,cTW,omegaTW,NuTW,thetaAmp,dphi];

    if Nx==32 && Nz==15, spanBase32=struct('S',S,'F',F2,'c',cSigned); end
    if Nx==48 && Nz==23, spanBase48=struct('S',S,'F',F2,'c',cSigned); end
    if Nx==64 && Nz==31, spanBase64=struct('S',S,'F',F2,'c',cSigned); end
    if Nx==80 && Nz==39, spanBase80=struct('S',S,'F',F2,'c',cSigned); end
    if saveSnapshots
        allHigh=struct('S',S,'F',F2,'tt',tt,'rr',rr,'snaps',{snaps},'snapTimes',snapTimes);
    end
end

T=array2table(rows,'VariableNames',{'Nx','Nz','dt','omega_standing','Nu_standing', ...
    'wave_speed_abs','omega_traveling','Nu_traveling','theta11_amp','phase_phi_minus_theta_abs'});
writetable(T,fullfile(outdir,'independent_grid_convergence_matlab.csv'));

% Galerkin reference appended separately for direct comparison.
Tref=table(NaN,NaN,NaN,1.339791,1.008716,0.400020,0.400020*sqrt(11.5155204), ...
    1.013882,0.586293,NaN,'VariableNames',T.Properties.VariableNames);
writetable([T;Tref],fullfile(outdir,'independent_grid_convergence_matlab_with_reference.csv'));

% High-grid transition series.
H=allHigh;
Tout=array2table([H.tt,H.rr],'VariableNames',{'t','odd_L2_norm','critical_fft_phase', ...
    'critical_fft_amplitude','Nu','phase_phi_over_theta'});
writetable(Tout,fullfile(outdir,'independent_transition_64x31_matlab.csv'));

% Co-moving invariance using the stored late snapshots.
[relAll,relTheta,cLate]=comovingResidual(H.S,H.snaps,H.snapTimes);
Tc=table(H.snapTimes(:),relAll(:),relTheta(:),'VariableNames', ...
    {'t','relative_L2_all_fields','relative_L2_theta'});
writetable(Tc,fullfile(outdir,'comoving_invariance_64x31_matlab.csv'));
fprintf('Late co-moving speed = %.9f; max relative L2 = %.3e\n',cLate,max(relAll));

% Compare against independently verified values; tolerances are intentionally
% broader than roundoff because MATLAB and SciPy FFT/DST paths need not be
% bitwise identical.
refFile=fullfile(outdir,'independent_grid_convergence.csv');
if isfile(refFile)
    V=readtable(refFile,'TextType','string');
    for k=1:min(height(T),height(V)-1)
        dc=100*abs(T.wave_speed_abs(k)-V.wave_speed_abs(k))/V.wave_speed_abs(k);
        dNu=100*abs(T.Nu_traveling(k)-V.Nu_traveling(k))/V.Nu_traveling(k);
        fprintf('%dx%d: speed diff %.4g%%, Nu diff %.4g%%\n',T.Nx(k),T.Nz(k),dc,dNu);
    end
else
    fprintf(['Optional independent reference CSV not found: %s\n' ...
             'Skipping cross-implementation percentage comparison; MATLAB outputs are still generated.\n'],refFile);
end

%% Independent time-step refinement at fixed 48 x 23 spatial resolution
% These runs support the manuscript statement that the remaining difference
% is dominated by spatial rather than temporal discretization.
dtList=[0.0100 0.0075 0.0050]; dtRows=zeros(numel(dtList),6);
for id=1:numel(dtList)
    dt=dtList(id); S=makeSolver(48,23,dt,1.0);
    F=initialField(S,0.70,0);
    [F,ts,rec]=advanceGrid(S,F,80,true,0.04,false);
    [X,Z]=ndgrid(S.x,S.z);
    for j=1:3
        F(:,:,j)=F(:,:,j)+1e-4*real(S.q(j))*sin(S.a*X).*sin(pi*Z);
    end
    [F2,tt,rr]=advanceGrid(S,F,160,false,0.04,false);
    late=tt>100; ph=-unwrap(rr(:,2)); pp=polyfit(tt(late),ph(late),1);
    dtRows(id,:)=[dt,abs(pp(1)/S.a),abs(pp(1)),mean(rr(late,4)), ...
        critAmplitude(S,F2,1),rr(end,1)];
end
Tdt=array2table(dtRows,'VariableNames',{'dt','wave_speed_abs','omega_traveling', ...
    'Nu_traveling','theta11_amp','odd_norm_final'});
writetable(Tdt,fullfile(outdir,'independent_timestep_refinement_48x23_matlab.csv'));

%% Second forcing level: deltaRa = 0.5 on the independent grid
S=makeSolver(48,23,0.0075,0.5); F=initialField(S,0.52,0);
[F,ts,rec]=advanceGrid(S,F,90,true,0.04,false);
[X,Z]=ndgrid(S.x,S.z);
for j=1:3
    F(:,:,j)=F(:,:,j)+1e-4*real(S.q(j))*sin(S.a*X).*sin(pi*Z);
end
[F2,tt,rr]=advanceGrid(S,F,220,false,0.04,false);
late=tt>160; ph=-unwrap(rr(:,2)); pp=polyfit(tt(late),ph(late),1);
row2=[0.5,48,23,0.0075,abs(pp(1)/S.a),abs(pp(1)),mean(rr(late,4)), ...
    critAmplitude(S,F2,1),rr(end,1)];
T2=array2table(row2,'VariableNames',{'deltaRa','Nx','Nz','dt','wave_speed_abs', ...
    'omega_traveling','Nu_traveling','theta11_amp','odd_norm_final'});
writetable(T2,fullfile(outdir,'second_parameter_dRa_0p5_matlab.csv'));

%% Limited 3-D spanwise secondary stability of the traveling state
% Linear disturbances are proportional to exp(i*k_y*y) in the frame moving
% with the measured traveling-wave speed. Positive sigma3D denotes growth.
% The refined 64x31 and 80x39 calculations are concentrated near the weak
% long-wave neutral crossing, where spatial resolution is most important.
kyRatio=[0.10 0.15 0.18 0.19 0.192 0.194 0.1955 0.1956 0.1957 0.196 0.198 0.20 0.25 0.50 1.00 1.50 2.00]';
spanRows=nan(numel(kyRatio),5);
set64=[0.15 0.18 0.19 0.192 0.194 0.196 0.198 0.20];
set80=[0.15 0.192 0.194 0.1955 0.1956 0.1957 0.196 0.20];
for ik=1:numel(kyRatio)
    rr=kyRatio(ik);
    s48=spanwiseGrowth(spanBase48.S,spanBase48.F,spanBase48.c,rr*spanBase48.S.a,0.00375,80,2);
    s32=NaN; s64=NaN; s80=NaN;
    if ismembertol(rr,[0.10 0.15 0.18 0.20 0.25 0.50 1.00],1e-12)
        s32=spanwiseGrowth(spanBase32.S,spanBase32.F,spanBase32.c,rr*spanBase32.S.a,0.005,120,2);
    end
    if ismembertol(rr,set64,1e-12)
        s64=spanwiseGrowth(spanBase64.S,spanBase64.F,spanBase64.c,rr*spanBase64.S.a,0.00375,80,2);
    end
    if ismembertol(rr,set80,1e-12)
        s80=spanwiseGrowth(spanBase80.S,spanBase80.F,spanBase80.c,rr*spanBase80.S.a,0.00375,80,2);
    end
    spanRows(ik,:)=[rr,s32,s48,s64,s80];
    fprintf('3-D k_y/a=%g: sigma32=% .6e sigma48=% .6e sigma64=% .6e sigma80=% .6e\n', ...
        rr,s32,s48,s64,s80);
end
T3=array2table(spanRows,'VariableNames',{'ky_over_a','growth_exponent_32x15', ...
    'growth_exponent_48x23','growth_exponent_64x31','growth_exponent_80x39'});
writetable(T3,fullfile(outdir,'limited_3D_spanwise_stability_matlab.csv'));
writetable(T3,fullfile(root,'limited_3D_spanwise_stability_matlab.csv'));

% Near-neutral time-step refinement on the finest 80 x 39 base state.
dt3=[0.00375 0.00250]; ky3=0.1956; dt3rows=zeros(numel(dt3),3);
for id=1:numel(dt3)
    sig=spanwiseGrowth(spanBase80.S,spanBase80.F,spanBase80.c,ky3*spanBase80.S.a,dt3(id),80,2);
    dt3rows(id,:)=[ky3,dt3(id),sig];
end
T3dt=array2table(dt3rows,'VariableNames',{'ky_over_a','dt','growth_exponent_80x39'});
writetable(T3dt,fullfile(outdir,'limited_3D_timestep_refinement_matlab.csv'));
writetable(T3dt,fullfile(root,'limited_3D_timestep_refinement_matlab.csv'));

%% ---------------- local functions --------------------------------------
function S=makeSolver(Nx,Nz,dt,dRa)
    S.Nx=Nx; S.Nz=Nz; S.dt=dt;
    S.H=10; S.J=0.5; S.M=2; S.A=0.13164; S.gamma=0.3/0.7;
    S.Rac=52.4889464007; S.Ra=S.Rac+dRa; S.R=sqrt(S.Ra);
    S.a=sqrt(11.5155206); S.Lx=2*pi/S.a;
    S.dx=S.Lx/Nx; S.dz=1/(Nz+1);
    S.x=(0:Nx-1)'*S.dx; S.z=(1:Nz)*S.dz;
    S.kx=(2*pi/S.Lx)*[0:Nx/2-1 -Nx/2:-1]';
    n=1:Nz; S.lz=-4*sin(n*pi/(2*(Nz+1))).^2/S.dz^2;
    S.lap=S.lz+(-S.kx.^2);  % implicit expansion Nx x Nz
    S.alpha=(1/S.M)*(1-S.J/S.M);
    j=(1:Nz)'; nn=1:Nz;
    S.Sz=sqrt(2/(Nz+1))*sin(j*nn*pi/(Nz+1));

    S.inv1=zeros(Nx,Nz,3,3); S.inv2=zeros(Nx,Nz,3,3);
    I=eye(3);
    for ix=1:Nx
        for iz=1:Nz
            lam=S.lap(ix,iz);
            buoy=S.Ra*S.kx(ix)^2/(-lam);
            L=[buoy+lam-S.H, S.H, 0; ...
               S.gamma*S.H/S.A, ((S.J/S.M)*lam-S.gamma*S.H)/S.A, S.alpha*lam/S.A; ...
               0, 1, -1/S.M];
            S.inv1(ix,iz,:,:)=inv(I-dt*L);
            S.inv2(ix,iz,:,:)=inv(3*I-2*dt*L);
        end
    end

    Lam=S.a^2+pi^2;
    Lc=[S.Rac*S.a^2/Lam-Lam-S.H, S.H, 0; ...
        S.gamma*S.H/S.A, -((S.J/S.M)*Lam+S.gamma*S.H)/S.A, -S.alpha*Lam/S.A; ...
        0,1,-1/S.M];
    [V,D]=eig(complex(Lc)); ev=diag(D); ids=find(imag(ev)>0);
    [~,jj]=min(abs(real(ev(ids)))); q=V(:,ids(jj)); q=q/norm(q);
    q=q*exp(-1i*angle(q(1))); S.q=q;
end

function F=initialField(S,amp,odd)
    [X,Z]=ndgrid(S.x,S.z); F=zeros(S.Nx,S.Nz,3);
    for j=1:3
        F(:,:,j)=amp*real(S.q(j))*cos(S.a*X).*sin(pi*Z);
        if odd~=0
            F(:,:,j)=F(:,:,j)+odd*real(S.q(j))*sin(S.a*X).*sin(pi*Z);
        end
    end
end

function C=tr(S,F)
    C=zeros(size(F));
    for j=1:3, C(:,:,j)=fft(F(:,:,j),[],1)*S.Sz; end
end

function F=itr(S,C)
    F=zeros(size(C));
    for j=1:3, F(:,:,j)=real(ifft(C(:,:,j)*S.Sz',[],1)); end
end

function Y=apply3(A,X)
    Y=zeros(size(X));
    for i=1:3
        tmp=zeros(size(X,1),size(X,2));
        for j=1:3, tmp=tmp+squeeze(A(:,:,i,j)).*X(:,:,j); end
        Y(:,:,i)=tmp;
    end
end

function [Nh,u,w]=Nhat(S,F)
    th=F(:,:,1); thh=fft(th,[],1)*S.Sz;
    wh=(-S.R*(S.kx.^2)./S.lap).*thh; wh(1,:)=0;
    w=real(ifft(wh*S.Sz',[],1));
    wz=zeros(size(w));
    wz(:,2:end-1)=(w(:,3:end)-w(:,1:end-2))/(2*S.dz);
    wz(:,1)=w(:,2)/(2*S.dz); wz(:,end)=-w(:,end-1)/(2*S.dz);
    wzh=fft(wz,[],1); uh=zeros(size(wzh));
    for ix=2:S.Nx, uh(ix,:)=1i*wzh(ix,:)/S.kx(ix); end
    u=real(ifft(uh,[],1));
    thx=real(ifft((1i*S.kx).*fft(th,[],1),[],1));
    thz=zeros(size(th));
    thz(:,2:end-1)=(th(:,3:end)-th(:,1:end-2))/(2*S.dz);
    thz(:,1)=th(:,2)/(2*S.dz); thz(:,end)=-th(:,end-1)/(2*S.dz);
    N=zeros(size(F)); N(:,:,1)=-(u.*thx+w.*thz); Nh=tr(S,N);
end

function F=symmetrize(F)
    Nx=size(F,1); idx=[1 Nx:-1:2]; F=0.5*(F+F(idx,:,:));
end

function [Fout,tout,rec,snaps,snapTimes]=advanceGrid(S,F,T,enforceSym,saveDt,saveSnapshots)
    dt=S.dt; nsteps=round(T/dt); saveEvery=max(1,round(saveDt/dt));
    Y0=tr(S,F); [N0,u,w]=Nhat(S,F);
    Y1=apply3(S.inv1,Y0+dt*N0); F1=itr(S,Y1);
    if enforceSym, F1=symmetrize(F1); Y1=tr(S,F1); end
    [N1,u,w]=Nhat(S,F1);
    maxSave=floor(nsteps/saveEvery)+2; tout=zeros(maxSave,1); rec=zeros(maxSave,5); ns=0;
    snaps={}; snapTimes=[]; snapEvery=max(1,round(0.5/dt));
    [ns,tout,rec]=recordState(S,F,0,tout,rec,ns);
    [ns,tout,rec]=recordState(S,F1,dt,tout,rec,ns);
    Ym=Y0; Yn=Y1; Nm=N0; Nn=N1; Fn=F1;
    for k=1:nsteps-1
        Yp=apply3(S.inv2,4*Yn-Ym+2*dt*(2*Nn-Nm)); Fp=itr(S,Yp);
        if enforceSym, Fp=symmetrize(Fp); Yp=tr(S,Fp); end
        [Np,u,w]=Nhat(S,Fp);
        if mod(k+1,saveEvery)==0
            [ns,tout,rec]=recordState(S,Fp,(k+1)*dt,tout,rec,ns,w);
        end
        if saveSnapshots && mod(k+1,snapEvery)==0 && (k+1)*dt>T-15
            snaps{end+1}=Fp; %#ok<AGROW>
            snapTimes(end+1,1)=(k+1)*dt; %#ok<AGROW>
        end
        Ym=Yn; Yn=Yp; Nm=Nn; Nn=Np; Fn=Fp;
    end
    Fout=Fn; tout=tout(1:ns); rec=rec(1:ns,:);
end

function [ns,tout,rec]=recordState(S,F,t,tout,rec,ns,w)
    if nargin<7, [~,~,w]=Nhat(S,F); end
    C=tr(S,F); ct=C(2,1,1); cp=C(2,1,2);
    ns=ns+1; tout(ns)=t;
    rec(ns,:)=[oddNorm(S,F),angle(ct),abs(ct),NuGrid(S,F,w),angle(cp/ct)];
end

function n=oddNorm(S,F)
    idx=[1 S.Nx:-1:2]; O=0.5*(F-F(idx,:,:));
    n=sqrt(S.dx*S.dz*sum(O(:).^2));
end

function Nu=NuGrid(S,F,w)
    avg=S.dz*mean(sum(w.*F(:,:,1),2));
    Nu=1+avg/(S.R*(1+1/S.gamma));
end

function [C,Sv]=critProjection(Svct,F,field)
    [X,Z]=ndgrid(Svct.x,Svct.z); f=F(:,:,field);
    C=4*Svct.dz/Svct.Nx*sum(f.*cos(Svct.a*X).*sin(pi*Z),'all');
    Sv=4*Svct.dz/Svct.Nx*sum(f.*sin(Svct.a*X).*sin(pi*Z),'all');
end

function amp=critAmplitude(S,F,field)
    [c,s]=critProjection(S,F,field); amp=hypot(c,s);
end

function omega=fitOscillationFrequency(t,y,wlo,whi)
    t=t(:); y=y(:); y=y-mean(y);
    obj=@(w) sinusoidResidual(w,t,y);
    omega=fminbnd(obj,wlo,whi,optimset('TolX',1e-10,'Display','off'));
end

function r=sinusoidResidual(w,t,y)
    A=[cos(w*t),sin(w*t),ones(size(t))];
    c=A\y; e=y-A*c; r=e'*e;
end

function y=wrapPi(x), y=mod(x+pi,2*pi)-pi; end

function [relAll,relTheta,c]=comovingResidual(S,snaps,t)
    n=numel(snaps); ph=zeros(n,1);
    for k=1:n
        C=tr(S,snaps{k}); ph(k)=-angle(C(2,1,1));
    end
    ph=unwrap(ph); pp=polyfit(t,ph,1); c=pp(1)/S.a;
    refphi=ph(1); aligned=cell(n,1);
    for k=1:n
        delta=(ph(k)-refphi)/S.a; Fh=fft(snaps{k},[],1);
        phase=exp(1i*S.kx*delta);
        for j=1:3, Fh(:,:,j)=phase.*Fh(:,:,j); end
        aligned{k}=real(ifft(Fh,[],1));
    end
    ref=aligned{1}; base=sqrt(S.dx*S.dz*sum(ref(:).^2));
    th=ref(:,:,1); baseTh=sqrt(S.dx*S.dz*sum(th(:).^2));
    relAll=zeros(n,1); relTheta=zeros(n,1);
    for k=1:n
        d=aligned{k}-ref; relAll(k)=sqrt(S.dx*S.dz*sum(d(:).^2))/base;
        dth=aligned{k}(:,:,1)-th; relTheta(k)=sqrt(S.dx*S.dz*sum(dth(:).^2))/baseTh;
    end
end


function sigma=spanwiseGrowth(S,Fbase,c,ky,dt,T,seed)
    % Leading exponent of a single spanwise Fourier class exp(i*k_y*y).
    [~,ub,wb]=Nhat(S,Fbase);
    thb=Fbase(:,:,1);
    thbx=ifft((1i*S.kx).*fft(thb,[],1),[],1);
    thbz=zeros(size(thb));
    thbz(:,2:end-1)=(thb(:,3:end)-thb(:,1:end-2))/(2*S.dz);
    thbz(:,1)=thb(:,2)/(2*S.dz); thbz(:,end)=-thb(:,end-1)/(2*S.dz);

    lap3=S.lz+(-S.kx.^2)-ky^2;
    kh2=S.kx.^2+ky^2;
    inv1=zeros(S.Nx,S.Nz,3,3); inv2=zeros(S.Nx,S.Nz,3,3); I=eye(3);
    for ix=1:S.Nx
        for iz=1:S.Nz
            lam=lap3(ix,iz); buoy=S.Ra*kh2(ix)/(-lam); adv=1i*c*S.kx(ix);
            L=[buoy+lam-S.H+adv, S.H, 0; ...
               S.gamma*S.H/S.A, ((S.J/S.M)*lam-S.gamma*S.H)/S.A+adv, S.alpha*lam/S.A; ...
               0,1,-1/S.M+adv];
            inv1(ix,iz,:,:)=inv(I-dt*L); inv2(ix,iz,:,:)=inv(3*I-2*dt*L);
        end
    end

    rng(seed,'twister'); G=randn(S.Nx,S.Nz,3)+1i*randn(S.Nx,S.Nz,3);
    G=G/spanNorm(S,G); Y0=trComplex(S,G); N0=spanwiseExplicit(S,G,ub,wb,thbx,thbz,ky);
    Y1=apply3(inv1,Y0+dt*N0); G1=itrComplex(S,Y1); N1=spanwiseExplicit(S,G1,ub,wb,thbx,thbz,ky);
    Ym=Y0; Yn=Y1; Nm=N0; Nn=N1;
    nsteps=round(T/dt); ren=max(1,round(1/dt)); acc=0; nt=floor(nsteps/ren); tt=zeros(nt,1); ll=zeros(nt,1); ii=0;
    for k=1:nsteps-1
        Yp=apply3(inv2,4*Yn-Ym+2*dt*(2*Nn-Nm)); Gp=itrComplex(S,Yp);
        Np=spanwiseExplicit(S,Gp,ub,wb,thbx,thbz,ky);
        Ym=Yn; Yn=Yp; Nm=Nn; Nn=Np;
        if mod(k+1,ren)==0
            nrm=spanNorm(S,Gp); acc=acc+log(nrm);
            Yn=Yn/nrm; Ym=Ym/nrm; Nn=Nn/nrm; Nm=Nm/nrm;
            ii=ii+1; tt(ii)=(k+1)*dt; ll(ii)=acc;
        end
    end
    tt=tt(1:ii); ll=ll(1:ii); use=tt>=0.4*T; pp=polyfit(tt(use),ll(use),1); sigma=pp(1);
end

function Nhat3=spanwiseExplicit(S,G,ub,wb,thbx,thbz,ky)
    th=G(:,:,1);
    thx=ifft((1i*S.kx).*fft(th,[],1),[],1);
    thz=zeros(size(th),'like',th);
    thz(:,2:end-1)=(th(:,3:end)-th(:,1:end-2))/(2*S.dz);
    thz(:,1)=th(:,2)/(2*S.dz); thz(:,end)=-th(:,end-1)/(2*S.dz);
    [du,~,dw]=spanwiseVelocity(S,th,ky);
    N=zeros(size(G),'like',G); N(:,:,1)=-(ub.*thx+wb.*thz+du.*thbx+dw.*thbz);
    Nhat3=trComplex(S,N);
end

function [u,v,w]=spanwiseVelocity(S,th,ky)
    thh=fft(th,[],1)*S.Sz;
    lap3=S.lz+(-S.kx.^2)-ky^2; kh2=S.kx.^2+ky^2;
    wh=(-S.R*kh2./lap3).*thh; w=ifft(wh*S.Sz',[],1);
    wz=zeros(size(w),'like',w);
    wz(:,2:end-1)=(w(:,3:end)-w(:,1:end-2))/(2*S.dz);
    wz(:,1)=w(:,2)/(2*S.dz); wz(:,end)=-w(:,end-1)/(2*S.dz);
    wzh=fft(wz,[],1); uh=(1i*S.kx).*wzh./kh2; vh=(1i*ky)*wzh./kh2;
    u=ifft(uh,[],1); v=ifft(vh,[],1);
end

function C=trComplex(S,F)
    C=zeros(size(F),'like',F);
    for j=1:3, C(:,:,j)=fft(F(:,:,j),[],1)*S.Sz; end
end

function F=itrComplex(S,C)
    F=zeros(size(C),'like',C);
    for j=1:3, F(:,:,j)=ifft(C(:,:,j)*S.Sz',[],1); end
end

function n=spanNorm(S,F)
    n=sqrt(S.dx*S.dz*sum(abs(F(:)).^2));
end
