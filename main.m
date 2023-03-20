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

reflectScrIons = true;
reflectAccelIons = false;
plotReflection = false;
trackParticles = true;

dt = 1.0035e-9;
Ntime = 2e-6/dt;
% noNewParticles = 100;

N_r = 301;
N_z = 751;
N = N_r*N_z;

Lz = 0.005;
Lr = 0.002;

Tplasma = 5;
RhoPlasma = 1.2e+17;
s_IonCharge = 1.60217657*1e-19;
s_IonMass = 131.293*1.66053892*1e-27;
mp_num = 1000;     % macro particle number for ions
mp_mass = s_IonMass*mp_num;
mp_charge = s_IonCharge*mp_num;
vBohm = sqrt(s_IonCharge * Tplasma / s_IonMass);
vRadial = sqrt(s_IonCharge * 0.04 / s_IonMass);
noNewParticles = uint8(dt*RhoPlasma*vBohm*(Lr^2 * 0.5*cosd(30)*sind(30))/mp_num);

dz = Lz/(N_z-1);
dr = Lr/(N_r-1);

Vplasma = 2266;
Vscreen = 2241;
Vaccel = -400;
Vplume = 0;

rAccel = 0.4e-3;
rScreen = 0.8e-3;
wAccel = 0.8e-3;
wScreen = 0.4e-3;
z0Screen = 1.0e-3;
dzAccelScreen = 1.0e-3;

ScreenBeginNodeRadial = rScreen/dr + 1;
ScreenBeginNodeAxial = z0Screen/dz + 1;
ScreenEndNodeRadial = N_r;
ScreenEndNodeAxial = (z0Screen+wScreen)/dz + 1;
patchScreenX = dz*[ScreenBeginNodeAxial ScreenEndNodeAxial ScreenEndNodeAxial ScreenBeginNodeAxial];
patchScreenY = dr*[ScreenBeginNodeRadial ScreenBeginNodeRadial ScreenEndNodeRadial ScreenEndNodeRadial];

AccelBeginNodeRadial = rAccel/dr + 1;
AccelBeginNodeAxial = (z0Screen+wScreen+dzAccelScreen)/dz + 1;
AccelEndNodeRadial = N_r;
AccelEndNodeAxial = (z0Screen+wScreen+dzAccelScreen+wAccel)/dz + 1;
patchAccelX = dz*[AccelBeginNodeAxial AccelEndNodeAxial AccelEndNodeAxial AccelBeginNodeAxial];
patchAccelY = dr*[AccelBeginNodeRadial AccelBeginNodeRadial AccelEndNodeRadial AccelEndNodeRadial];

x_screen = ScreenBeginNodeAxial*dz;
x_accel = AccelBeginNodeAxial*dz;

z0 = 1e-6;

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
V = zeros(N,1);
potential = zeros(N_r, N_z);

tic
i = 1;
% Configure Coefficient Matrix
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
t_ConfigureLoop = toc
% Construct sparse matrix A from the index and value vectors
A = sparse(idx,idy,val);
tic
% Since majority of A is 0, converting A to sparse matrix will be much more
% efficient. Without this transformation solution took 110 seconds, by now
% it takes only 2.5 seconds
[Potential,Ex,Ey] = Solve(A,b,N_r,N_z);
t_PotentialSolve = toc

particles(1) = particle([z0,Lr*0.3], [vBohm,0], mp_charge, mp_mass, trackParticles);
particles(2) = particle([z0,Lr*0.5], [vBohm,0], mp_charge, mp_mass, trackParticles);
particles(3) = particle([z0,Lr*0.8], [vBohm,0], mp_charge, mp_mass, trackParticles);

tic
for i = 1:Ntime
    % Create new particles at the inlet with random velocities according to
    % Maxvellian velocity distribution with mean of vBohm
%     if mod(i,50) == 0
    if i == 1
        % Create number of particle objects
        new(1:noNewParticles) = particle;

        % Assign each new object some random initial values
        for j = 1:length(new)
            velx = vBohm;
            vely = vRadial*randn;
            while (vely < -4*vRadial || vely > 4*vRadial)
                vely = vRadial*randn;
            end
            new(j) = particle([z0 rand*Lr], [velx vely], mp_charge, mp_mass, trackParticles);
        end
        particles = [particles new];
    end
    hold off 

    % Initialize array to store position and trajectories
    X = zeros(1,length(particles));
    Y = zeros(1,length(particles));
    if trackParticles
        trajX = zeros(i,length(particles));
        trajY = zeros(i,length(particles));
    end
    
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

        % Check if the particle left the simulation domain
        if particles(k).pos(2) <= 0
            particles(k) = particles(k).ReflectBottom();
        elseif particles(k).pos(2) >= Lr
            particles(k) = particles(k).ReflectTop(Lr);
        elseif particles(k).pos(1) >= Lz
            particles(k) = [];
            continue
        elseif particles(k).pos(1) < 0
            particles(k) = [];
            continue
        end
        
        % Check if the particle hit any grid
        if (posX > z0Screen && posX < (z0Screen+wScreen) && posY > rScreen && posY < Lr)
            if reflectScrIons
                particles(k) = particles(k).ReflectScreen(x_screen);
            else 
                particles(k) = [];
                continue
            end
        end       
        if (posX > (z0Screen+wScreen+dzAccelScreen) && posX < (z0Screen+wScreen+dzAccelScreen+wAccel) && posY > rAccel && posY < Lr)
            if reflectAccelIons
                particles(k) = particles(k).ReflectAccel(x_accel);
            else
                particles(k) = [];
                continue
            end
        end

        % Caltulate the index of the the bounding nodes
        if (particles(k).pos(1) >= 0) && (particles(k).pos(2) >= 0)
            idx_lower = round(particles(k).pos(1)/dz - 0.5) +1;
            idy_lower = round(particles(k).pos(2)/dr - 0.5) +1;
            idx_upper = ceil(particles(k).pos(1)/dz) +1;
            idy_upper = ceil(particles(k).pos(2)/dr) +1;
        else
            disp('SOMETHING WRONG')
        end
        % Calculate the areas with the bounding nodes for density and
        % electric field weighting
        A1 = (posX - (idx_lower-1)*dz)*(posY - (idy_lower-1)*dr); % lower-lower area
        A2 = ((idx_upper-1)*dz - posX)*(posY - (idy_lower-1)*dr); % upper-lower area
        A3 = (posX - (idx_lower-1)*dz)*((idy_upper-1)*dr - posY); % lower-upper area 
        A4 = ((idx_upper-1)*dz - posX)*((idy_upper-1)*dr - posY); % upper-upper area 

        interpEx = (A4*Ex(idy_lower,idx_lower) + A3*Ex(idy_lower,idx_upper)...
                  + A2*Ex(idy_upper,idx_lower) + A1*Ex(idy_upper,idx_upper))/(dr*dz);
        interpEy = (A4*Ey(idy_lower,idx_lower) + A3*Ey(idy_lower,idx_upper)...
                  + A2*Ey(idy_upper,idx_lower) + A1*Ey(idy_upper,idx_upper))/(dr*dz);
        % Calculate the gradient and move the particle
        gradient = mp_num*[interpEx interpEy];
        particles(k) = particles(k).MoveInField(gradient,dt);

        % Store position and trajectory data in matrices to plot later on.
        % Plotting by a single plot function is faster
        X(k) = particles(k).pos(1);
        Y(k) = particles(k).pos(2);        
        if trackParticles
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
        k = k + 1;
    end 
    numParticle = length(particles);  % number of particles in the simulation domain
    cmap = jet(length(particles));    % create color array for each particle
    % Plotting each particle in a loop and holding on is rather a slow
    % process, if we store these data and plot it in a single turn it is
    % much faster
    scatter(X(1:numParticle),Y(1:numParticle),20,cmap,'filled'), hold on
    % Patch() function is substantially faster than plot(), however I
    % couldn't manage to color each line with a different RGB value
    if trackParticles
        patch(trajX(:,1:numParticle), trajY(:,1:numParticle),'green')  % here 'green' is needed to satisfy number of input arguments but not represented in figure
    end
    patch(patchScreenX,patchScreenY,'red')    
    patch(patchAccelX,patchAccelY,'blue')    
    if plotReflection
        scatter(X(1:numParticle),-Y(1:numParticle),20,cmap,'filled')
        patch(patchAccelX,-patchAccelY,'blue')
        patch(patchScreenX,-patchScreenY,'red')
        if trackParticles
            patch(trajX(:,1:numParticle), -trajY(:,1:numParticle),'green')
        end
        ylim([-Lr Lr])
    else
        ylim([0 Lr])
    end
    title(sprintf('%d MacroParticles In the System -- Completed: %.2f %%', numParticle,100*i/Ntime))
    xlim([0 Lz])  
    
    drawnow
end
toc

figure(2)
surf(1:N_z, 1:N_r, Potential,'LineStyle','--')
hold on
surf(1:N_z, linspace(1,-N_r,N_r), Potential,'LineStyle','--')



