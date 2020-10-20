%This requires a Comsol Multiphysics license with LiveLink for Matlab (no
%further Comsol toolboxes required)
%More details on code are given in sweep_density_Pe_params.m
k=10;
load(['parameters/p',num2str(k),'.mat'])
rho = 0.1;
    Pe_sweep = [logspace(2,3,100),1500, 2000, 2500, 3000, 3500, 4000];
for i = 1:size(Pe_sweep,2)
    i
        %t_max helps find the correct steady solution. Depends on parameters.
        t_max = 1e2;
        [~,b(i),c(i)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,Pe_sweep(i),rho,t_max,1);
        if c(i)<0
            disp('Rescaling t_max')
            t_max = t_max/10;
            [~,b(i),c(i)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,Pe_sweep(i),rho,t_max,1);
        end
        b(i)
        c(i)
end
%Output results
        save('restime/steady.mat')
        figure
        loglog(Pe_sweep,b)
        savefig('restime/steady.fig')

        