%This code generates data for steady 2D simulations used in Fig. 1
%This requires a Comsol Multiphysics 5.5 license with LiveLink for Matlab (no
%further Comsol toolboxes required)


%Write parameter files (can comment this out if they're already written)
%write_parameters;
N = length(dir('parameters/*.mat'));
%Sweep over parameter files
for k = 1:N
clearvars -except k
%Load parameter file
load(['parameters/p',num2str(k),'.mat'])
shear_rate_sweep = [100 500 1000 5000 10000]; %in 1/s
Pe_sweep = p.H * p.H * shear_rate_sweep/(500); %p.H in microns, D in microns^2/s. This is a Peclet number based on the size of the computational domain being the same as the height of the cell layer, with u=1 at the top.
rhoc = zeros(size(Pe_sweep));
%Write k to terminal to keep track
k
%To get the steady solution we perform a time-dependent simulation for
%t_max in non-dimensional time, and use it as an initial condition.
%t_max may need to be different for different parameter sets - it is
%calibrated automatically but only by reduction to save time - if nothing works then
%increase this initial value by a factor of 10
t_max = 1e4; %use t_max = 1e4 for typical parameters
%initialise variable for tracking any failed attempts
fail_counter = [];
%Sweep over Peclet numbers
for j = 1:size(Pe_sweep,2)
%Write j to terminal to keep track
j
p.Pe = Pe_sweep(j);
Pe_eff(j) = 4.062* p.Pe^(1/3) / (3^(5/3)  * p.L^(1/3)); %effective Peclet number (used to predict critical density)
%Define parameter groups
        Lambda = p.lambda * p.mu / (p.alpha * p.gamma);
        K = p.beta * (p.km+p.gamma) / (p.kp*p.gamma);
            %Prediction of critical density
            options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-12,'StepTolerance',1e-12);
            fun = @(x) (x)*tan(x)- Pe_eff(j);
            omega_guess = fsolve(fun,pi/4,options);
            rho_guess = K*(omega_guess^2 + p.kappa)/(p.r*Lambda);
            %order of magnitude of rho_guess (for sweeping over rho)
            n = floor(log(abs(rho_guess))./log(10));
%Sweep values of rho in log space (around rho_guess)
rho_sweep{j} = rho_guess / 10^(n) *logspace(n-1/10,n+1/10,100);
%Skip if predicted rhoc>1, since we are underestimating the minimum
%effective Peclet number so the numerical rhoc>1 too
            if rho_guess>1
                disp('Rhoc greater than 1. Skipping')
                continue
            end
%Resize domain if necessary (Peclet number must be scaled appropriately if so because we still
%apply u=1 at the top of the domain)
if p.L/Pe_sweep(j)>0.1
H_3 = 10;
p.Pe = 10*p.Pe;
else
H_3 = 1;
end
%Calibrate t_max to save time and avoid negative steady solutions
disp('Calibrating tmax')
[~,b{j}(100),c{j}(100)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,p.Pe,rho_sweep{j}(100),t_max,H_3);
    if c{j}(100)<0
        disp('Rescaling tmax')
        t_max = t_max/10;
        [~,b{j}(100),c{j}(100)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,p.Pe,rho_sweep{j}(100),t_max,H_3);
        if c{j}(100)<0
            disp('Rescaling tmax again')
            t_max = t_max/10;
            [~,b{j}(100),c{j}(100)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,p.Pe,rho_sweep{j}(100),t_max,H_3);
            if c{j}(100)<0
                disp('Rescaling tmax again again')
                t_max = t_max/10;
                [~,b{j}(100),c{j}(100)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,p.Pe,rho_sweep{j}(100),t_max,H_3);
                if c{j}(100)<0
                    disp('Rescaling tmax again again again')
                    t_max = t_max/10;
                    [~,b{j}(100),c{j}(100)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,p.Pe,rho_sweep{j}(100),t_max,H_3);
                    if c{j}(100)<0
                        disp('Did not find appropriate tmax: skipping')
                        continue
                    end
                end
            end
        end
    end
    %If not on upper branch, need to check for higher densities. Then
    %calibrate t_max again
    if b{j}(100)<1
        disp('Max density still on lower branch - trying for higher densities')
        rho_sweep{j} = rho_guess / 10^(n) *logspace(n,n+1/5,100);
        %Edit t_max if necessary
        disp('Calibrating tmax')
        [~,b{j}(100),c{j}(100)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,p.Pe,rho_sweep{j}(100),t_max,H_3);
            if c{j}(100)<0
             disp('Rescaling tmax')
                t_max = t_max/10;
             [~,b{j}(100),c{j}(100)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,p.Pe,rho_sweep{j}(100),t_max,H_3);
                if c{j}(100)<0
                 disp('Rescaling tmax again')
                    t_max = t_max/10;
                 [~,b{j}(100),c{j}(100)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,p.Pe,rho_sweep{j}(100),t_max,H_3);
                 if c{j}(100)<0
                     disp('Did not find appropriate tmax: skipping')
                    continue
            end
        end
    end
    end
    %Abandon if still on lower branch and return later (use higher densities)
    b{j}(100)
    c{j}(100)
    if b{j}(100)<1
        disp('Still on lower branch - return to check what is happening later')
        continue
    end
        
%Calculate rhoc by sweeping over rho
for i = 1:size(rho_sweep{j},2)
    %Keep track of progress - comment out to take up less space in terminal
    if mod(i,10)==0
        i
    end
    %Output time of simulation (comment out tic and toc to save space in
    %terminal)
    %tic
    %Perform simulation
        [~,b{j}(i),c{j}(i),d{j}(i),x_num{j}{i}(1,:),Pe_num{j}{i}(1,:),Cmax{j}(i)] = solve_steady_comsol(p.q,p.r,p.lambda,p.mu,p.kp,p.km,p.alpha,p.beta,p.gamma,p.kappa,p.L,p.Pe,rho_sweep{j}(i),t_max,H_3);
    %toc
    %Check the first solution is not already on the upper branch
    if i==1
        if b{j}(i)>1
            disp('WARNING: first solution already on upper branch. Use a lower value of the density.')
            break
        end
    end
    if c{j}(i)<0
        disp('WARNING: Negative solution. Use a better initial condition.')
        break
    end
end
    %Skip to next value if loop over rho was not completed
    if size(b{j},2)<size(rho_sweep{j},2)
        omegac(j) = 0;
        fail_counter = [fail_counter j];
        continue
    end
        %Output figure (saved to results folder later)
        figure
        loglog(rho_sweep{j},b{j})
        hold on
        loglog(rho_sweep{j},c{j})
        %Now find the critical density
        %This requires curve fitting toolbox. If it is not available, use
        %numerical gradient, commented out below (much less accurate)
        %bxx = gradient(gradient(b)./gradient(rho_sweep))./gradient(rho_sweep);
        s = spline(log(rho_sweep{j}),log(b{j}));
        x = linspace(min(rho_sweep{j}),max(rho_sweep{j}),1000);
        z = ppval(s,log(x));
        hold on;
        loglog(x,exp(z),'--')
        p_der=fnder(s,2);
        bxx{j} = ppval(p_der,log(x));
        [~,idx] = max(bxx{j});
        rhoc(j) = x(idx);
        plot([rhoc(j) rhoc(j)],[10^(-3) 10^(2)],'--')
        savefig(['res/p',num2str(k),'j',num2str(j),'.fig'])
        close
        omegac(j) = sqrt((rhoc(j) * p.r * Lambda - p.kappa * K) /(K));
        Pec(j) = omegac(j)*tan(omegac(j));

end

save(['res/p',num2str(k),'.mat'])
end       