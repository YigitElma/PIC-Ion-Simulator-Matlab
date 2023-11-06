% -------------------------------------------------------------------------
%     Yigit Gunsur ELMACIOGLU
%     15.03.2023
% 
% Ion Thruster Grid Region Simulator
% 
% -------------------------------------------------------------------------
clc
clear
close all

% Simulation Settings -----------------------------------------------------
bool_reflectScrIons = true;      % true: reflect ions from the screen grid as ions
bool_reflectAccelIons = false;   % true: reflect ions from the accel grid as ions
bool_plotReflection = false;     % true: plots the mirror image of the solution
bool_trackParticles = true;      % true: stores and plots particle trajectories
bool_solvePoisson = true;       % true: solve Poisson's eqn false: Laplace
bool_saveAnimation = false;      % true: records the plots to save .avi file (might use too much memory)
freq_Poisson = 1;     % frequency for solving Poisson eqn
freq_Plot = 1000;       % frequency for plotting recent solution
freq_IonCreate = 1;   % frequency to create new ion particles

% The percentage of kinetic energy preserved after reflection from a grid (between 0-1)
percentEnergyAfterReflection = 0.1;
% Survival probability of an ion to be reflected as an ion from the screen grid (between 0-1)
survivalProb = 0.1; 
% Number of real particles represented by a macro-particle
mp_num = 100;     

% Mesh Properties ---------------------------------------------------------
N_r = 301;
N_z = 751;
N  = N_r*N_z;
Lz = 0.005;
Lr = 0.002;
dz = Lz/(N_z-1);
dr = Lr/(N_r-1);

% Electric Potential Values in Volt ---------------------------------------
Vplasma = 2266;
Vscreen = 2241;
Vaccel  = -400;
Vplume  = 0;

Tplasma = 5;       % Temperature of the plasma in eV
Tneutral = 0.04;   % Temperature of the neutral part in eV
RhoPlasma = 1.2e+17;     % ion number density of the plasma
s_IonCharge = 1.60217657*1e-19;        % single ion electric charge
s_IonMass = 131.293*1.66053892*1e-27;  % single Xenon atom mass
mp_mass = s_IonMass*mp_num;      % mass of a macro-particle
mp_charge = s_IonCharge*mp_num;  % electric charge of a macro-particle
vBohm = sqrt(s_IonCharge * Tplasma / s_IonMass);     % Bohm Velocity
vRadial = sqrt(s_IonCharge * Tneutral / s_IonMass);  % Neutral gas thermal velocity

ev2joule = 1.6e-19;     % conversion constant from eV to Joules
% Choose dt such that when particle have max kinetic energy,
% it doesn't move a distance greater than mesh size in a single time step
dt = dz/sqrt((Tplasma + Vplasma - Vaccel)*ev2joule*2/s_IonMass); 
Ntime = 0.4e-6/dt;    % Number of simulation steps for time advancing
% Number of incoming ions (as macro-particles) are calculated assuming the inlet is a 30-90-60
% triangle, and ions enter with Bohm velocity from a plasma which has ion
% density of RhoPlasma ( # = flux*dt = (density)*(velocity)*(area)*dt )
noNewParticles = uint8(dt*RhoPlasma*vBohm*(Lr^2 * 0.5*cosd(30)*sind(30))/mp_num);

% Geometric Properties of the Simulation Domain in meter ------------------
rAccel   = 0.4e-3;
rScreen  = 0.8e-3;
wAccel   = 0.8e-3;
wScreen  = 0.4e-3;
z0Screen = 1.0e-4;
dzAccelScreen = 1.0e-3;

% Node indices for the grids ----------------------------------------------
ScreenBeginNodeRadial = rScreen/dr + 1;
ScreenBeginNodeAxial  = z0Screen/dz + 1;
ScreenEndNodeRadial   = N_r;
ScreenEndNodeAxial    = (z0Screen+wScreen)/dz + 1;

AccelBeginNodeRadial = rAccel/dr + 1;
AccelBeginNodeAxial  = (z0Screen+wScreen+dzAccelScreen)/dz + 1;
AccelEndNodeRadial   = N_r;
AccelEndNodeAxial    = (z0Screen+wScreen+dzAccelScreen+wAccel)/dz + 1;

% For plotting the screen and accel grid on the plots ---------------------
patchScreenX = dz*[ScreenBeginNodeAxial ScreenEndNodeAxial ScreenEndNodeAxial ScreenBeginNodeAxial];
patchScreenY = dr*[ScreenBeginNodeRadial ScreenBeginNodeRadial ScreenEndNodeRadial ScreenEndNodeRadial];
patchAccelX  = dz*[AccelBeginNodeAxial AccelEndNodeAxial AccelEndNodeAxial AccelBeginNodeAxial];
patchAccelY  = dr*[AccelBeginNodeRadial AccelBeginNodeRadial AccelEndNodeRadial AccelEndNodeRadial];

% Axial position of the screen and accel grid upstream face in meter ------
z_screen = z0Screen;
z_accel  = z0Screen+wScreen+dzAccelScreen;
z0 = 1e-8;   % to prevent weighting errors, give non-zero z0, but very small value

% Constructing A as a sparse matrix reduces the memory usage drastically,
% to be able to solve larger N, you should use sparse from the
% beginning, otherwise 32Gb of memory may not be enough
% A = sparse(N, N);
% sparse(N,N) is good but the configuration might take longer since A is
% sparse, if you use spalloc(N,N,x) it will preallocate x number of nonzero
% elements inside A sparse matrix and the configureation will take less for
% larger N values
% A = spalloc(N,N,Nsparse);
% To accelerate further, calculate number of non-zero elements ahead
NsparseNorm = 5*(N_r-2)*(N_z-2) + 2*(N_r+N_z-2);
NsparseGrid = (ScreenEndNodeRadial-ScreenBeginNodeRadial)*(ScreenEndNodeAxial-ScreenBeginNodeAxial+1)...
              +(AccelEndNodeRadial-AccelBeginNodeRadial)*(AccelEndNodeAxial-AccelBeginNodeAxial+1);
Nsparse = NsparseNorm - 4*NsparseGrid + 2*(N_z-2) - (ScreenEndNodeAxial-ScreenBeginNodeAxial+1) - (AccelEndNodeAxial-AccelBeginNodeAxial+1);
% Store i index, j index and the value of the elements in seperate vectors,
% we will use these to construct sparse A later
idx = zeros(Nsparse,1);
idy = zeros(Nsparse,1);
val = ones(Nsparse,1);
b = zeros(N,1);

i = 1;
% Configure Coefficient Matrix --------------------------------------------
for z = 1:N_z
    for r = 1:N_r
        index = (r-1)*N_z + z;
        if z == 1
            % A(index, index) = 1;
            idx(i) = index;
            idy(i) = index;
            val(i) = 1;
            b(index) = Vplasma;
        elseif z == N_z
            % A(index, index) = 1;
            idx(i) = index;
            idy(i) = index;
            val(i) = 1;
            b(index) = Vplume;
        elseif (z >= ScreenBeginNodeAxial && z <= ScreenEndNodeAxial && r >= ScreenBeginNodeRadial && r <= ScreenEndNodeRadial)
            % A(index, index) = 1;
            idx(i) = index;
            idy(i) = index;
            val(i) = 1;
            b(index) = Vscreen;
        elseif (z >= AccelBeginNodeAxial && z <= AccelEndNodeAxial && r >= AccelBeginNodeRadial && r <= AccelEndNodeRadial)
            % A(index, index) = 1;
            idx(i) = index;
            idy(i) = index;
            val(i) = 1;
            b(index) = Vaccel;
        elseif r == 1
            % A(index, index) = 1;
            idx(i) = index;
            idy(i) = index;
            val(i) = 1;
            i = i+1;
            % A(index, index+N_z) = -1;
            idx(i) = index;
            idy(i) = index+N_z;
            val(i) = -1;
            b(index) = 0;
        elseif r == N_r
            % A(index, index) = 1;            
            idx(i) = index;
            idy(i) = index;
            val(i) = 1;
            i = i+1;
            % A(index, index-N_z) = -1;
            idx(i) = index;
            idy(i) = index-N_z;
            val(i) = -1;
            b(index) = 0;
        else
            % A(index, index) = -2*( 1/dz^2 + 1/dr^2 );
            idx(i) = index;
            idy(i) = index;
            val(i) = -2*( 1/dz^2 + 1/dr^2 );
            i = i+1;
            % A(index, index-1) = 1/dz^2 ;
            idx(i) = index;
            idy(i) = index-1;
            val(i) = 1/dz^2;
            i = i+1;
            % A(index, index+1) = 1/dz^2 ;
            idx(i) = index;
            idy(i) = index+1;
            val(i) = 1/dz^2;
            i = i+1;
            % A(index, index+N_z) = 1/dr^2 + 1/2/dr;
            idx(i) = index;
            idy(i) = index+N_z;
            val(i) = 1/dr^2 + 1/2/dr;
            i = i+1;
            % A(index, index-N_z) = 1/dr^2 - 1/2/dr;
            idx(i) = index;
            idy(i) = index-N_z;
            val(i) = 1/dr^2 - 1/2/dr;
        end
        i = i+1;
    end
end

% Construct sparse matrix A from the index and value vectors
A = sparse(idx,idy,val);
% Since majority of A is 0, converting A to sparse matrix will be much more
% efficient. Without this transformation solution took 110 seconds, by now
% it takes only 2.5 seconds
[Potential,Ez,Er] = Solve(A,b,N_r,N_z,dz,dr);

% Create some initial particles (number 3 is random here)
particles(1) = particle([z0 rand*Lr], [vBohm,0], mp_charge, mp_mass, bool_trackParticles);
particles(2) = particle([z0 rand*Lr], [vBohm,0], mp_charge, mp_mass, bool_trackParticles);
particles(3) = particle([z0 rand*Lr], [vBohm,0], mp_charge, mp_mass, bool_trackParticles);
  
t = 1;
fig = figure;
fig.WindowState = 'maximized';

for i = 1:Ntime
    % Create new particles at the inlet with random velocities according to
    % Maxvellian velocity distribution with mean of vBohm
    if mod(i,freq_IonCreate) == 0
%     if i == 1      % use this for fast trajectory check
        % Create number of particle objects
        new(1:noNewParticles) = particle;
        % Assign each new object some random initial values
        for j = 1:length(new)
            velx = vBohm;
            vely = vRadial*randn;
            while (vely < -4*vRadial || vely > 4*vRadial)
                vely = vRadial*randn;
            end
            new(j) = particle([z0 rand*Lr], [velx vely], mp_charge, mp_mass, bool_trackParticles);
        end
        particles = [particles new];
    end
    hold off 

    % Initialize array to store position and trajectories
    X = zeros(1,length(particles));
    Y = zeros(1,length(particles));
    if bool_trackParticles
        trajX = zeros(i,length(particles));
        trajY = zeros(i,length(particles));
    end

    RhoIons = zeros(N,1);
    
    % Since some particles are deleted during the loop, I used 'while'
    % because it allows the upper boundary of the loop to be changed
    % during the loop. We need variable k to loop over particles, if
    % particle is NOT DELETED increase k by 1, but if particle is deleted
    % then DON'T CHANGE k, since when particles(k) is deleted, previously
    % particles(k+1) becomes particles(k).
    k = 1;
    while k <= length(particles)    
        posX = particles(k).pos(1);
        posY = particles(k).pos(2);

        % Caltulate the index of the the bounding nodes
        if (particles(k).pos(1) >= 0) && (particles(k).pos(2) >= 0)
            idx_lower = round(posX/dz - 0.5) +1;
            idy_lower = round(posY/dr - 0.5) +1;
            idx_upper = ceil(posX/dz) +1;
            idy_upper = ceil(posY/dr) +1;
        else
            disp('SOMETHING WRONG')
        end
        % Calculate the areas with the bounding nodes for density and
        % electric field weighting
        A1 = (posX - (idx_lower-1)*dz)*(posY - (idy_lower-1)*dr); % lower-lower area
        A2 = ((idx_upper-1)*dz - posX)*(posY - (idy_lower-1)*dr); % upper-lower area
        A3 = (posX - (idx_lower-1)*dz)*((idy_upper-1)*dr - posY); % lower-upper area 
        A4 = ((idx_upper-1)*dz - posX)*((idy_upper-1)*dr - posY); % upper-upper area 

        interpEz = (A4*Ez(idy_lower,idx_lower) + A3*Ez(idy_lower,idx_upper)...
                  + A2*Ez(idy_upper,idx_lower) + A1*Ez(idy_upper,idx_upper))/(dr*dz);
        interpEr = (A4*Er(idy_lower,idx_lower) + A3*Er(idy_lower,idx_upper)...
                  + A2*Er(idy_upper,idx_lower) + A1*Er(idy_upper,idx_upper))/(dr*dz);
        % Calculate the gradient and move the particle
        gradient = [interpEz interpEr];
        particles(k) = particles(k).MoveInField(gradient,dt);

        % Update the position variables
        posX = particles(k).pos(1);
        posY = particles(k).pos(2);

        % Check if the particle left the simulation domain
        if posY <= 0
            particles(k) = particles(k).ReflectBottom();
        elseif posY >= Lr
            particles(k) = particles(k).ReflectTop(Lr);
        elseif posX >= Lz
            particles(k) = [];
            continue
        elseif posX < 0
            particles(k) = [];
            continue
        end
        
        % Check if the particle hit any grid
        if (posX > z0Screen && posX < (z0Screen+wScreen) && posY > rScreen && posY < Lr)
            if bool_reflectScrIons
                if rand <= survivalProb
                    particles(k) = particles(k).ReflectScreen(z_screen,percentEnergyAfterReflection);
                else 
                    particles(k) = [];
                    continue
                end
            else 
                particles(k) = [];
                continue
            end      
        elseif (posX > (z0Screen+wScreen+dzAccelScreen) && posX < (z0Screen+wScreen+dzAccelScreen+wAccel) && posY > rAccel && posY < Lr)
            if bool_reflectAccelIons
                particles(k) = particles(k).ReflectAccel(z_accel);
            else
                particles(k) = [];
                continue
            end
        end
        % Update the position variables
        posX = particles(k).pos(1);
        posY = particles(k).pos(2);

        % Caltulate the index of the the bounding nodes
        if (particles(k).pos(1) >= 0) && (particles(k).pos(2) >= 0)
            idx_lower = round(posX/dz - 0.5) +1;
            idy_lower = round(posY/dr - 0.5) +1;
            idx_upper = ceil(posX/dz) +1;
            idy_upper = ceil(posY/dr) +1;
        else
            disp('SOMETHING WRONG')
        end
        % Calculate the areas with the bounding nodes for density and
        % electric field weighting
        A1 = (posX - (idx_lower-1)*dz)*(posY - (idy_lower-1)*dr); % lower-lower area
        A2 = ((idx_upper-1)*dz - posX)*(posY - (idy_lower-1)*dr); % upper-lower area
        A3 = (posX - (idx_lower-1)*dz)*((idy_upper-1)*dr - posY); % lower-upper area 
        A4 = ((idx_upper-1)*dz - posX)*((idy_upper-1)*dr - posY); % upper-upper area 

        % Charge densities at nodes for Poisson Solution
        if bool_solvePoisson
            index = (idy_lower-1)*N_z + idx_lower;
            RhoIons(index) = RhoIons(index) + mp_charge/(dr*dz^2)*A4/(dr*dz);  % I just assumed cell volume to be dr*dz^2
            RhoIons(index+1) = RhoIons(index+1) + mp_charge/(dr*dz^2)*A3/(dr*dz);
            RhoIons(index+N_z) = RhoIons(index+N_z) + mp_charge/(dr*dz^2)*A2/(dr*dz);
            RhoIons(index+N_z+1) = RhoIons(index+N_z+1) + mp_charge/(dr*dz^2)*A1/(dr*dz);
        end

        if mod(i,freq_Plot) == 0
            % Store position and trajectory data in matrices to plot later on.
            % Plotting by a single plot function is faster
            X(k) = particles(k).pos(1);
            Y(k) = particles(k).pos(2);
            if bool_trackParticles
                num = length(particles(k).trajectory(1,:));
                if num == i
                    trajX(:,k) = particles(k).trajectory(1,:);
                    trajY(:,k) = particles(k).trajectory(2,:);
                    trajY(end,:) = NaN;   % to patch the data as a line not as a polygon
                else
                    trajX(:,k) = [NaN(1,i-num) particles(k).trajectory(1,:)];
                    trajY(:,k) = [NaN(1,i-num) particles(k).trajectory(2,:)];
                end
            end
        end
        k = k + 1;
    end 

    if mod(i,freq_Plot) == 0
        numParticle = length(particles);  % number of particles in the simulation domain
        cmap = jet(length(particles));    % create color array for each particle
        % Plotting each particle in a loop and holding on is rather a slow
        % process, if we store these data and plot it in a single turn it is
        % much faster
        scatter(X(1:numParticle),Y(1:numParticle),20,cmap,'filled'), hold on
        % Patch() function is substantially faster than plot(), however I
        % couldn't manage to color each line with a different RGB value
        if bool_trackParticles
            patch(trajX(:,1:numParticle), trajY(:,1:numParticle),'green')  % here 'green' is needed to satisfy number of input arguments but not represented in figure
        end
        patch(patchScreenX,patchScreenY,'red')    
        patch(patchAccelX,patchAccelY,'blue')    
        if bool_plotReflection
            patch(patchAccelX,-patchAccelY,'blue')
            patch(patchScreenX,-patchScreenY,'red')
            if bool_trackParticles
                patch(trajX(:,1:numParticle), -trajY(:,1:numParticle),'green')
            end
            scatter(X(1:numParticle),-Y(1:numParticle),20,cmap,'filled')
            ylim([-Lr Lr])
        else
            ylim([0 Lr])
        end
        title(sprintf('%d MacroParticles In the System -- Completed: %.2f %%', numParticle,100*i/Ntime))
        xlim([0 Lz])        
        if bool_saveAnimation
            anim(t) = getframe;
            t = t+1; 
        end       
        drawnow
    end

%     % Solve Poisson's eqn and update potential and gradient fields
%     if bool_solvePoisson && mod(i,freq_Poisson) == 0
%         bPoisson = RhoIons + b;
%         [Potential,Ez,Er] = Solve(A,bPoisson,N_r,N_z,dz,dr);
%     end
end

% Create .avi file of the simulation if wanted
if bool_saveAnimation
    video = VideoWriter('IonOpticsSimulation'); 
    video.FrameRate = 100;
    open(video);
    writeVideo(video, anim);
    close(video);
end

% Plot the electric potential distribution  as a valley
figure(2)
surf((0:N_z-1)*dz, (0:N_r-1)*dr, Potential,'LineStyle','--')
hold on
surf((0:N_z-1)*dz, linspace(0,(-N_r+1)*dr,N_r), Potential,'LineStyle','--')



