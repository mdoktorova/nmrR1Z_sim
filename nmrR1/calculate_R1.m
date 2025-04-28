%-------------------------------------------------------------------------
% Method for calculating CH bond R1Z relaxation rates from MD simulations
% 
% Protocol based on Doktorova, Khelashvili, Ashkar and Brown, 2022
% https://doi.org/10.1016/j.bpj.2022.12.007
%
% Milka Doktorova, November 2022
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% calculate_R1.m
%
% Uses the ACF functions (mangleacfJi_X.dat files) of all carbons defined
% in ch_vectors.txt to calculate the corresponding R1 relaxation rates
%
% NOTE: uses parallel computing toolbox in MATLAB. If not available,
% replace 'parfor' with 'for'.
%-------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Input files
%
% ch_vectors.txt
%           first line has the name of the lipid (resname in VMD's
%           notation); each subsequent line has three atom names C H1 H2
%           corresponding to a carbon whose two CH bonds are to be
%           analyzed. For double bonds where a carbon has a single hydrogen
%           atom, the name of the hydrogen atom should appear twice
%
% mangleacfJ0_X.txt, mangleacfJ1_X.txt, mangleacfJ2_X.txt
%           X is the carbon name (e.g. C22, C23, C34, C35, etc); the file
%           has the ACFs G0, G1 and G2 of the CH bonds at this carbon
%           (averaged over the two carbon's hydrogen atoms specified in
%           ch_vectors.txt); it has NL columns and NF rows where NL is the
%           number of lipids that have this carbon and NF is number of lags
%           of the ACF (usually total number of simulation frames divided
%           by 2)
%
% -------------------------------------------------------------------------
% Output files (in results_dtt directory)
% 
% Jomega_J0_X.dat; Jomega_J1_X.dat; Jomega_J2_X.dat 
%           calculated spectral densities for the full range of sampled
%           frequencies; file has two columns corresponding to the Larmor
%           frequency in MHz and the corresponding spectral density
%
% spectral_J0_X.dat; spectral_J1_CX.dat; spectral_J2_CX.dat
%           a b c coefficients of best power-law fits to the respective 
%           spectral densities (a*x^b+c)
%
% R1_X.dat
%           R1Z relaxation rate calculated for carbon X in 1/s
%
% coeffs.dat
%           b coefficients from best power-law fits to J0, J1 and J2; file
%           has a row for each carbon specified in ch_vectors.txt with 3 
%           entries being the b coefficients for the fits to J0, J1 and J2 
%           spectral densities of that carbon
%
% dtt_allG.dat
%           smallest timesteps found for resampling the G0, G1 and G2 ACFs;
%           file has a row for each carbon specified in ch_vectors.txt with
%           3 entries being the smallest dt values for G0, G1 and G2 ACFs 
%           of that caborn. Note that all calculations proceed by  
%           resampling all 3 ACFs with the smallest timestep found for G0
%
% rsq.dat
%           R-square goodness of fit for the best fits to the J0, J1 and J2
%           spectral densities; file has 1 line with 3 entries (3 R-square 
%           parameters) for the fits to J0, J1 and J2 for each carbon 
%           specified in ch_vectors.txt, in that order
%
%-------------------------------------------------------------------------

% -----------------------------------------------------
% USER-specified input parameters
% -----------------------------------------------------

% output frequency of trajectory frames in microseconds
dt = 40*(10^(-6));

% Larmor frequency in MHz
LF = 76.8;

% lowest frequency to be sampled as f/T MHz where T is total sampled time
% in microseconds (note: usually T is total time of the ACF equal to total 
% time of the simulation divided by 2); default is 1/T (so f=1)
f = 1;

% write/save all files: yes (default, 1) or no (2)
wflag = 1;

% -----------------------------------------------------
% END OF USER-specified parameters
% -----------------------------------------------------

% read atom names from file
data = readcell("ch_vectors.txt");
carbons = [];
nc = length(data);
for i = 2:nc
    carbons = [carbons string(data{i})];
end

% load first ACF for reference
mangleacfJ0 = load(strcat('mangleacfJ0_',carbons(1),'.dat'));

% calculate total simulation time in microseconds
Nfrall = length(mangleacfJ0);
Nfr = ceil(Nfrall/2);
T = (Nfr-1)*dt;

% define the frequencies to be sampled
% from (f/T) to 1/(2*dt) MHz
maxF = 1/(2*dt);
r1nu = (f/T):(1/T):maxF;
omega = 2*pi*r1nu;

% initial time range in microseconds
t=[0:dt:T-dt]';

R1all = zeros(length(carbons),length(r1nu));
R1allc = zeros(length(carbons),1);
scd = []; rsq = [];
scddiff = []; scddiffFit = [];
coeffs = []; dttallG = [];

for ci=1:length(carbons)

    const = 4.2785*(10^(-2)); % 1/us^2
    
    % load ACF    
    mangleacfJ0 = readtable(strcat('mangleacfJ0_',carbons(ci),'.dat'));
    mangleacfJ1 = readtable(strcat('mangleacfJ1_',carbons(ci),'.dat'));
    mangleacfJ2 = readtable(strcat('mangleacfJ2_',carbons(ci),'.dat'));

    % set variance to the value of ACF at t=0
    v0 = abs(mangleacfJ0.Var1(Nfr));
    v1 = abs(mangleacfJ1.Var1(Nfr));
    v2 = abs(mangleacfJ2.Var1(Nfr));

    % get the ACF without element at t=0 and fit to power function
    mangleacfJ0 = abs(mangleacfJ0.Var1(Nfr+1:end));
    mangleacfJ1 = abs(mangleacfJ1.Var1(Nfr+1:end));
    mangleacfJ2 = abs(mangleacfJ2.Var1(Nfr+1:end));
    
    tfit = [dt:dt:T];
    [bestfit gof] = fit(tfit',mangleacfJ0,'power2');
    cval = coeffvalues(bestfit);
    a0c = cval(1); b0c = cval(2); c0c = cval(3);
    cfJ0 = @(x) a0c*(x.^b0c)+c0c;
    [bestfit gof] = fit(tfit',mangleacfJ1,'power2');
    cval = coeffvalues(bestfit);
    a1c = cval(1); b1c = cval(2); c1c = cval(3);
    cfJ1 = @(x) a1c*(x.^b1c)+c1c;
    [bestfit gof] = fit(tfit',mangleacfJ2,'power2');
    cval = coeffvalues(bestfit);
    a2c = cval(1); b2c = cval(2); c2c = cval(3);
    cfJ2 = @(x) a2c*(x.^b2c)+c2c;

    % find smallest dt at which the ACF at t=dt is less than the variance
    dtt_all = dt/400:dt/400:dt;
    dtt_comp = find(cfJ0(dtt_all)<=v0);
    dtt0 = dtt_all(dtt_comp(1));

    dtt_all = dt/400:dt/400:dt;
    dtt_comp = find(cfJ1(dtt_all)<=v1);
    dtt1 = dtt_all(dtt_comp(1));

    dtt_all = dt/400:dt/400:dt;
    dtt_comp = find(cfJ2(dtt_all)<=v2);
    dtt2 = dtt_all(dtt_comp(1));

    dttallG = [dttallG;dtt0 dtt1 dtt2];

    % define new (resampled) time range in microseconds
    dtt = dtt0; 
    t = dtt:dtt:T;

    % get the resampled ACFs
    mangleacfJ0 = cfJ0(t);
    mangleacfJ1 = cfJ1(t);
    mangleacfJ2 = cfJ2(t);

    % calculate spectral density from resampled ACF
    JomegaJ0 = zeros(length(r1nu),1);
    JomegaJ1 = zeros(length(r1nu),1);
    JomegaJ2 = zeros(length(r1nu),1);
    parfor i=1:length(r1nu)
        JomegaJ0(i) = 2*trapz(t,mangleacfJ0.*cos(omega(i)*t)) + v0*dtt;
        JomegaJ1(i) = 2*trapz(t,mangleacfJ1.*cos(omega(i)*t)) + v1*dtt;
        JomegaJ2(i) = 2*trapz(t,mangleacfJ2.*cos(omega(i)*t)) + v2*dtt;
    end

    % fit spectral density to a power function
    [bestfit gof] = fit(r1nu',JomegaJ0,'power2');
    cval = coeffvalues(bestfit);
    rsq = [rsq gof.rsquare];
    a0 = cval(1); b0 = cval(2); c0 = cval(3);
    if wflag == 1
        dlmwrite(strcat('results_dtt/spectral_J0_',carbons(ci),'.dat'),[a0 b0 c0],' ');
    end
    sdJ0 = @(x) a0*(x.^b0)+c0;
    
    [bestfit gof] = fit(r1nu',JomegaJ1,'power2');
    cval = coeffvalues(bestfit);
    rsq = [rsq gof.rsquare];
    a1 = cval(1); b1 = cval(2); c1 = cval(3);
    if wflag == 1
        dlmwrite(strcat('results_dtt/spectral_J1_',carbons(ci),'.dat'),[a1 b1 c1],' ');
    end
    sdJ1 = @(x) a1*(x.^b1)+c1;
    
    [bestfit gof] = fit(r1nu',JomegaJ2,'power2');
    cval = coeffvalues(bestfit);
    rsq = [rsq gof.rsquare];
    a2 = cval(1); b2 = cval(2); c2 = cval(3);
    if wflag == 1
        dlmwrite(strcat('results_dtt/spectral_J2_',carbons(ci),'.dat'),[a2 b2 c2],' ');
    end
    sdJ2 = @(x) a2*(x.^b2)+c2;
    
    % use fit to get carbon's relaxation rate at the specified Larmor Frequency
    R1allc(ci) = (const*((sdJ0(LF)+4*sdJ0(2*LF)) + 2*(sdJ1(LF)+4*sdJ1(2*LF)) + 2*(sdJ2(LF)+4*sdJ2(2*LF)))) * (10^6); % 1/s

    disp(carbons(ci));

    % write relaxation rate and spectral densities to files
    if wflag == 1
        dlmwrite(strcat('results_dtt/R1_',carbons(ci),'.dat'),R1allc(ci),' ');
    	dlmwrite(strcat('results_dtt/Jomega_J0_',carbons(ci),'.dat'),[r1nu' JomegaJ0],' ');
	    dlmwrite(strcat('results_dtt/Jomega_J1_',carbons(ci),'.dat'),[r1nu' JomegaJ1],' ');
	    dlmwrite(strcat('results_dtt/Jomega_J2_',carbons(ci),'.dat'),[r1nu' JomegaJ2],' ');
    end
    
    % save b coefficients of spectral densities J0, J1 and J2
    coeffs = [coeffs;b0 b1 b2];
    
end

if wflag == 1
    dlmwrite('results_dtt/coeffs.dat',coeffs,' ');
    dlmwrite('results_dtt/rsq.dat',rsq,' ');
    dlmwrite('results_dtt/dtt_allG.dat',dttallG,' ');
end

