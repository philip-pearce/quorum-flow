function [ p ] = parameters_comsol
%Dimensional parameters (subscript d is used)
Hd = 1; %microns
%(Note in our computations we take the cell population height to be 1 (i/e. it is scaled by Hd).
%However in these parameter files we set p.H = Hd, so that it is clear what height
%is being used (because this changes all the other parameters)
Dd = 500; %microns^2/s
Q0d = 5; %nM
rd = 0.01; %nM/s
kpd = 0.0001; %/nM/s
kmd = 0.01; %/s
qd = 1; %nM/s
lambdad = 0.01; %/s
mud = 1000; %/s
alphad = 0.001; %/s
betad = 0.01; %/s
gammad = 0.001; %/s
kappad = 0.0001; %/s

        % H is the domain height.
        p.H = 1;
        %base rate of R production
        p.r = Hd^2 * rd /(Dd * Q0d);
        %forward and backward binding of R and C, respectively
        p.kp = Hd^2 * kpd * Q0d /Dd;
        p.km = Hd^2 * kmd / Dd;
            
        % q is the base production rate of QSM
        p.q =  Hd^2 /(Dd * Q0d);
        
        % lambda is the coefficient of QSM production due to I
        p.lambda = Hd^2 * lambdad / Dd;
        %mu is the coefficient of I production due to C
        p.mu = Hd^2 * mud / Dd;
        
        % alpha, beta, and gamma are the degradation rates of I, R, and C,
        % respectively.
        p.alpha = Hd^2 * alphad / Dd;        
        p.beta = Hd^2 * betad / Dd;
        p.gamma =Hd^2 * gammad / Dd;
        %degradation of QS molecules
        p.kappa = Hd^2 * kappad / Dd;
    

end

