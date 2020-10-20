delete parameters/*.mat
delete parameters/*.txt
H_sweep = linspace(1,10,10);
for i =1:10
        [ p ] = parameters;
        H = H_sweep(i);
        %Multiply kinetic parameters by H^2
        p = structfun(@(x) H^2*x,p,'UniformOutput',false);
        %Correct H
        p.H = p.H/H;
%(Note in our computations we use the cell population height as Hd, so H=1 always.
%However in these parameter files we set p.H = Hd (the dimensional value), so that
%it is clear what height
%is being used (because this changes all the other parameters)
        %Aspect ratio 0.5
        p.L=2;
    writetable(struct2table(p),['parameters/p',num2str(i),'.txt'])
    save(['parameters/p',num2str(i),'.mat'],'p')
end