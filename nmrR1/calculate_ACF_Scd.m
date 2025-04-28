%-------------------------------------------------------------------------
% Method for calculating CH bond R1Z relaxation rates from MD simulations
% 
% Protocol based on Doktorova, Khelashvili, Ashkar and Brown, 2022
% https://doi.org/10.1016/j.bpj.2022.12.007
%
% Milka Doktorova, November 2022
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% calculate_ACF_Scd.m
%
% Uses the CH vectors (angle_X.dat files) of all carbons defined in 
% ch_vectors.txt to calculate the corresponding autocorrelation functions
% (ACF) and CH bond order parameters (Scd)
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
%           atom, the name of the same hydrogen atom should appear twice
%
% angle_X.txt
%           X is the carbon name (e.g. C22, C23, C34, C35, etc); the file
%           has NL x 6 columns and N rows where NL is the number of
%           lipids that have this carbon and N is the number of frames in
%           the trajectory. The 6 entries for each lipid x1 y1 z1 x2 y2 z2
%           are the directions of the two CH vectors at this carbon atom
%           obtained by subtracting the carbon coordinates from the
%           hydrogen coordinates
%
% -------------------------------------------------------------------------
% Output files
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
% scd_X.dat
%           the order parameter of carbon X
%

% -----------------------------------------------------
% no USER-specified input parameters
% -----------------------------------------------------

% read names of lipid and its carbbon and hydrogen atoms
data = readcell("ch_vectors.txt");
carbons = [];
nc = length(data);
for i = 2:nc
    carbons = [carbons string(data{i})];
end

for ci=1:length(carbons)
    disp('Carbon:');
    disp(ci);
    data = load(strcat('angle_',carbons(ci),'.txt'));

    % get the directions of the two CH vectors for all lipids    
    Nlip = length(data(1,:))/6;
    xa = data(:,1:Nlip);
    ya = data(:,1*Nlip+1:2*Nlip);
    za = data(:,2*Nlip+1:3*Nlip);

    xb = data(:,3*Nlip+1:4*Nlip);
    yb = data(:,4*Nlip+1:5*Nlip);
    zb = data(:,5*Nlip+1:6*Nlip);

    % initialize the arrays for the Wigner rotation matrix elements
    Nframes = length(xa(:,1));
    j0a = zeros(Nframes,Nlip);
    j1a = zeros(Nframes,Nlip);
    j2a = zeros(Nframes,Nlip);
    j0b=j0a; j1b=j1a; j2b=j2a;

    % define the bilayer normals for the top and bottom leaflets
    bntop = [0 0 1];
    bnbot = [0 0 -1];

    % calculate the Wigner rotation matrix elements for each lipid in every
    % frame
    parfor i=1:Nframes
    	for j=1:Nlip
            if j<=Nlip/2
                bn = bntop;
            else
                bn = bnbot;
            end
            vec = [xa(i,j) ya(i,j) za(i,j)];
            s1 = (norm(cross(vec,bn))/(norm(vec)*norm(bn)));
            c1 = (dot(vec,bn)/(norm(vec)*norm(bn)));
    	    v1 = cross(vec,bn); v2 = [1 0 0];
    	    gamma = atan2(norm(cross(v1,v2)),dot(v1,v2));
            j0a(i,j) = 1.5*c1^2-0.5;
            j1a(i,j) = sqrt(3/2)*s1*c1*exp(-gamma*sqrt(-1));
            j2a(i,j) = sqrt(3/8)*s1^2*exp(-2*gamma*sqrt(-1));

            vec = [xb(i,j) yb(i,j) zb(i,j)];
            s2 = (norm(cross(vec,bn))/(norm(vec)*norm(bn)));
            c2 = (dot(vec,bn)/(norm(vec)*norm(bn)));
    	    v1 = cross(vec,bn); v2 = [1 0 0];
    	    gamma = atan2(norm(cross(v1,v2)),dot(v1,v2));
            j0b(i,j) = 1.5*c2^2-0.5;
            j1b(i,j) = sqrt(3/2)*s2*c2*exp(-gamma*sqrt(-1));
            j2b(i,j) = sqrt(3/8)*s2^2*exp(-2*gamma*sqrt(-1));
        end
    end

    % calculate the ACF for each Wigner rotation matrix element
    parfor i=1:Nlip
        disp(i)
        angleacfaJ0(:,i) = xcorr(j0a(:,i),floor(Nframes/2),'unbiased')-mean(j0a(:,i))^2;
    	angleacfbJ0(:,i) = xcorr(j0b(:,i),floor(Nframes/2),'unbiased')-mean(j0b(:,i))^2;

    	angleacfaJ1(:,i) = xcorr(j1a(:,i),floor(Nframes/2),'unbiased');
        angleacfbJ1(:,i) = xcorr(j1b(:,i),floor(Nframes/2),'unbiased');

    	angleacfaJ2(:,i) = xcorr(j2a(:,i),floor(Nframes/2),'unbiased');
        angleacfbJ2(:,i) = xcorr(j2b(:,i),floor(Nframes/2),'unbiased');
    end

    % average the ACFs of the Wigner rotation matrix elements of the two CH
    % bonds for each carbon and write the results
    mangleacf = mean([mean(angleacfaJ0,2) mean(angleacfbJ0,2)],2);
    dlmwrite(strcat('mangleacfJ0_',carbons(ci),'.dat'),mangleacf,' ');
    
    mangleacf = mean([mean(angleacfaJ1,2) mean(angleacfbJ1,2)],2);
    dlmwrite(strcat('mangleacfJ1_',carbons(ci),'.dat'),mangleacf,' ');

    mangleacf = mean([mean(angleacfaJ2,2) mean(angleacfbJ2,2)],2);
    dlmwrite(strcat('mangleacfJ2_',carbons(ci),'.dat'),mangleacf,' ');

    % calculate and write the carbon's order parameter Scd
    scd = mean([mean(mean(j0a)); mean(mean(j0b))]);
    dlmwrite(strcat('scd_',carbons(ci),'.dat'),scd,' ');

end
