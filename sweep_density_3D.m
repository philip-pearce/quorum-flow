%This requires a Comsol Multiphysics license with LiveLink for Matlab (no
%further Comsol toolboxes required)
cd parameters
%We let H=5 microns
load('p5.mat')
cd ..
p.Pe = 50; %gammadot = 1000/s

%Sweep values of rho (choose manually)
rho_sweep = [[0.1:0.01:0.16], [0.161:0.001:0.18], [0.185:0.005:0.2]];
%Sweep over rho
for i = 1:size(rho_sweep,2)
    i
    tic
    %Perform steady solve (similar method to 2D case). Solve Comsol file to
    %disk for analysis later.
        [sol,b(i)] = solve_steady_comsol_3D(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.Pe,rho_sweep(i),10,2);
        rho_sweep(i)
        b
        sol.save(['res3D/sol_amp2_',num2str(i)])
    toc
end
        %Plot solution
        figure
        loglog(rho_sweep,b)


        