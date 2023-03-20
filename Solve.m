function [potential, Ex, Ey] = Solve(A,b,N_r,N_z)
    potential = zeros(N_r, N_z);
   
    % Solve linear system of equations    
    V = A\b;

    % Map 1D V vector to 2D
    for i = 1:N_r
        potential(i,:) = V( ((i-1)*N_z + 1):(i*N_z) );
    end
    
    % Calculate the gradient of potential
    [Ex,Ey] = gradient(potential);

    % Electric field is the negative of the gradient 
    Ex = -Ex;
    Ey = -Ey;
end