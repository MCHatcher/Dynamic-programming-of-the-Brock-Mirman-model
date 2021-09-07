% Value function Matlab code 
% Code originally written by Panayiotis Pourpourides, pourpouridesp@cf.ac.ac.uk
% This version by Michael Hatcher, m.c.hatcher@soton.ac.uk

clc; close; clear;

%----------------------------------------------------
% 1. Set number of grid points for k (nz) and z (nz)
%----------------------------------------------------
nk = 1000; % No. of grid points for k % a high number will make k appear continuous
nz = 2; % No. of states of technology % if you change this number you also need to change Z and PI and simulation below (see N states code)

%-----------------
% 2. Calibration
%-----------------
sigma = 2;
theta = 0.3;
delta = 0.1;
beta = 0.98;

%---------------------------
% 3. Steady state solution
%---------------------------
ky_ss = beta*theta/(1-beta*(1-delta));
y_ss = ky_ss^(theta/(1-theta));
k_ss = ky_ss * y_ss;
c_ss = y_ss - delta*k_ss;
z_ss = 1; % normalize the steady state for z to unity

%--------------------------------
% 4. Discretize the state space
%--------------------------------
%Two states
Z = [0.975 1.025];
PI = [0.6 0.4; 0.3 0.7];

Kstep = (1.2*k_ss-0.8*k_ss) / (nk-1); % increment for the k grid points
K = 0.8*k_ss:Kstep:1.2*k_ss; % the sequence of grid points starts from 0.8*k_ss and increases by Kstep until it reaches 1.2*k_ss

%------------------------------------
% 5. Set initial value function V(0)
%------------------------------------
u = c_ss^(1-sigma) / (1-sigma);
maxiter = 1000;
V = ones(nz,nk) * u / (1-beta); % This is the infinite sum of discounted utilities, when consumption remains always at the steady state.

%--------------------------------------
% 6. Iterate on V until convergence
%--------------------------------------
maxdiff = 1e-4; % = 0.0001
diffV = 1;
iter = 1;

while (iter <= maxiter) && (diffV > maxdiff) % as long as both statements hold the iteration continues
    
% Recall that V(k,z)=max{u(k,z,k')+ beta*EV(k',z')}
% Calculate expected V
%The first loop controls the rows of the matrix, the second loop controls the columns of the matrix and the third loop controls the sum in each cell.

EV = zeros(nz,nk);
for iz = 1:nz
    for ik_prime = 1:nk
        for iz_prime = 1:nz
    EV(iz,ik_prime) = EV(iz,ik_prime) + PI(iz,iz_prime) * V(iz_prime,ik_prime);
        end
    end
end

%--------------------
% 7. Compute new V
%--------------------
newV = zeros(nz,nk);
iDecK = zeros(nz,nk); % index to decisions in k grid
for iz = 1:nz
    for ik = 1:nk
    % Calculate c for each k'
    % c is a 1*nk vector
    
    c = Z(iz)*K(ik)^theta + (1-delta)*K(ik) - K;
    c = max(c,1e-8);
    
    % c must be non-negative. If there is a negative c replace it with a
    % very small number. Negative or small cs wouldn't be chosen anyway.
    
    % Compute the value function
    
    v = c.^(1-sigma)/(1-sigma) + beta*EV(iz,:); % This is a 1 x nk vector
    [newV(iz,ik), iDecK(iz,ik)] = max(v); % pick the highest V, and save it in newV in the corresponding location. 
    
    %iDecK(iz,ik) is the location of optimal k0 as a function of the location of the current state (z; k)
    end
end

%--------------------------------------------------
% 8. Compute the distance between V(n+1) and V(n)
%--------------------------------------------------
diffV = max(abs((newV(:)-V(:))./V(:)));
iter = iter + 1;
V = newV;

end

fprintf(1,'Found value function after %5.0d iterations.',iter-1)

%-----------------------------------
% 9. Plot the optimal rule for k'
%-----------------------------------
figure (1) % number the figures. If you omit this command, only the last figure will be displayed
plot(K,K(iDecK)')
title('Optimal decision rules for k-prime as a function of k. Different lines represent different z.')

%---------------------------------------------------------------
% 10. Simulate the model using the optimal rules for T periods
%---------------------------------------------------------------
T = 100;
simiZ = zeros(1,T);
simiK = zeros(1,T);
simiZ(1) = 1;
simiK(1) = 1;
rand('state',150);
for i = 2:T
    r = rand; % draw random number from uniform distribution on (0,1) interval
    if r<PI(simiZ(i-1),1), iz = 1; else iz=2;  
    end
ik = simiK(i-1); % timing: simK(t) = choice of k' in period t, i.e. k in period t+1
newiK = iDecK(iz,ik);
simiZ(i) = iz;
simiK(i) = newiK;
end

% There are 2 states for z, z1, z2. We used the uniform dist. on (0,1) to determine the realized state.  
 
%---------------------------------
% 11. Compute simulated y and c
%---------------------------------

simZ = Z(simiZ(2:end)); % a vector that contains the simulated value of z at each date
simKnext = K(simiK(2:end)); % a vector that contains the simulated grid point of k0 at each date
simK = K(simiK(1:end-1)); % a vector that contains the simulated grid point of k at each date
% y and c
simY = simZ.*simK.^theta;
simC = simY + (1-delta)*simK - simKnext;
%For distribution plot below
simZ = Z(simiZ(1:end));

%---------------------------------------
% 12. Plot simulated k, y and c
%---------------------------------------

figure (2)
subplot (221)
plot([simK'/k_ss])
title('Simulated k relative to its steady state')
subplot (222)
plot([simC'/c_ss])
title('Simulated c relative to its steady state')
subplot (223)
plot([simY'/y_ss])
title('Simulated y relative to its steady state')
subplot (224)
plot([(simY'-simC')/(y_ss - c_ss)])
title('Simulated s relative to its steady state') %saving = output - consumption as closed economy with no government
figure (3)
hist(simZ)
title('Productivity distribution. This will change with the PI entries and the number of states')