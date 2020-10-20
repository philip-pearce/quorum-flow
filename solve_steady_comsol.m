function [out,Qmax,Qmin,Qleft,x, Pe_eff,Cmax] = solve_steady_comsol(q,r,lambda,mu,kp,km,alpha,beta,gamma,kappa,L,Pe,rho,t_max,H_3)
%
% cleanest_2d.m
%
% Model exported on Sep 10 2020, 16:33 by COMSOL 5.5.0.359.

%Note there are some lines of Matlab code to change the mesh if H3 is not
%equal to 1. They are optimised for H3=10, and should be changed manually
% if any other value of H3 is used. If so, the Peclet number of the
% computation also needs to be scaled appropriately (we set u=1 at the top
% of the domain).

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.label('model.mph');

model.baseSystem('none');

%Define parameter values------------------------
model.param.set('q', num2str(q));
model.param.set('r', num2str(r));
model.param.set('lambda', num2str(lambda));
model.param.set('mu', num2str(mu));
model.param.set('kp', num2str(kp));
model.param.set('km', num2str(km));
model.param.set('alpha', num2str(alpha));
model.param.set('beta', num2str(beta));
model.param.set('gamma', num2str(gamma));
model.param.set('kappa', num2str(kappa));
model.param.set('Pe', num2str(Pe));
model.param.set('L', num2str(L));
%initial rho to use as IC - use guess for rhoc divided by 10, so on lower
%branch
model.param.set('rho', num2str(rho));
model.param.set('H3', num2str(H_3));
%------------------------------------------------

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').label('cells');
model.component('comp1').geom('geom1').feature('r2').set('pos', [0 0]);
model.component('comp1').geom('geom1').feature('r2').set('size', {'L' '1+H3'});
model.component('comp1').geom('geom1').create('pc1', 'ParametricCurve');
model.component('comp1').geom('geom1').feature('pc1').label('cell surface');
model.component('comp1').geom('geom1').feature('pc1').set('parmax','L');
model.component('comp1').geom('geom1').feature('pc1').set('coord', {'s' '1'});
model.component('comp1').geom('geom1').create('pard1', 'PartitionDomains');
model.component('comp1').geom('geom1').feature('pard1').set('partitionwith', 'edges');
model.component('comp1').geom('geom1').feature('pard1').selection('domain').set('r2(1)', 1);
model.component('comp1').geom('geom1').feature('pard1').selection('edge').set('r2(1)', 1);
model.component('comp1').geom('geom1').create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('r3').label('inflow');
model.component('comp1').geom('geom1').feature('r3').set('pos', [-5 1]);
model.component('comp1').geom('geom1').feature('r3').set('size', {'5' 'H3'});
model.component('comp1').geom('geom1').create('r4', 'Rectangle');
model.component('comp1').geom('geom1').feature('r4').label('outflow');
model.component('comp1').geom('geom1').feature('r4').set('pos', {'L' '1'});
model.component('comp1').geom('geom1').feature('r4').set('size',{'5' 'H3'});
model.component('comp1').geom('geom1').run;

model.variable.create('var1');
model.component('comp1').variable.create('var2');
model.component('comp1').variable('var2').set('u', '(y-1)/H3');
model.component('comp1').variable('var2').set('v', '0');
model.component('comp1').variable('var2').selection.geom('geom1', 2);
model.component('comp1').variable('var2').selection.set([1 3 4]);
model.component('comp1').variable.create('var3');
model.component('comp1').variable('var3').set('u', '0');
model.component('comp1').variable('var3').set('v', '0');
model.component('comp1').variable('var3').selection.geom('geom1', 2);
model.component('comp1').variable('var3').selection.set([2]);

model.component('comp1').physics.create('tds2', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds2').identifier('tds2');
model.component('comp1').physics('tds2').field('concentration').field('Q');
model.component('comp1').physics('tds2').field('concentration').component({'Q'});
model.component('comp1').physics('tds2').create('conc2', 'Concentration', 1);
model.component('comp1').physics('tds2').feature('conc2').selection.set([1]);
model.component('comp1').physics('tds2').create('reac1', 'Reactions', 2);
model.component('comp1').physics('tds2').feature('reac1').selection.set([2]);
model.component('comp1').physics('tds2').create('reac2', 'Reactions', 2);
model.component('comp1').physics('tds2').feature('reac2').selection.set([1 3 4]);
model.component('comp1').physics.create('tds', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds').field('concentration').field('I');
model.component('comp1').physics('tds').field('concentration').component({'I'});
model.component('comp1').physics('tds').selection.set([2]);
model.component('comp1').physics('tds').create('reac1', 'Reactions', 2);
model.component('comp1').physics('tds').feature('reac1').selection.set([2]);
model.component('comp1').physics.create('tds3', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds3').field('concentration').field('R');
model.component('comp1').physics('tds3').field('concentration').component({'R'});
model.component('comp1').physics('tds3').selection.set([2]);
model.component('comp1').physics('tds3').create('reac1', 'Reactions', 2);
model.component('comp1').physics('tds3').feature('reac1').selection.set([2]);
model.component('comp1').physics.create('tds4', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds4').field('concentration').field('C');
model.component('comp1').physics('tds4').field('concentration').component({'C'});
model.component('comp1').physics('tds4').selection.set([2]);
model.component('comp1').physics('tds4').create('reac1', 'Reactions', 2);
model.component('comp1').physics('tds4').feature('reac1').selection.set([2]);

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ref1', 'Refine');
model.component('comp1').mesh('mesh1').create('ref2', 'Refine');
if H_3<10
model.component('comp1').mesh('mesh1').create('ref3', 'Refine');

if H_3>1
    model.component('comp1').mesh('mesh1').create('ref4', 'Refine');
end
end

model.component('comp1').variable('var2').label('Flow field (fluid)');
model.component('comp1').variable('var3').label('Flow field (cells)');

model.component('comp1').view('view1').axis.set('xmin', -5.331976890563965);
model.component('comp1').view('view1').axis.set('xmax', 8.611043930053711);
model.component('comp1').view('view1').axis.set('ymin', -7.278668403625488);
model.component('comp1').view('view1').axis.set('ymax', 9.278668403625488);

model.component('comp1').physics('tds2').prop('MassConsistentStabilization').set('glim_mass', '0.1[mol/m^3]/tds2.helem');
model.component('comp1').physics('tds2').feature('cdm1').set('u', {'Pe*u'; 'Pe*v'; '0'});
model.component('comp1').physics('tds2').feature('cdm1').set('D_Q', [1; 0; 0; 0; 1; 0; 0; 0; 1]);
model.component('comp1').physics('tds2').feature('init1').set('initc', 'Q');
model.component('comp1').physics('tds2').feature('conc2').set('species', true);
model.component('comp1').physics('tds2').feature('reac1').set('R_Q', 'rho*(q+lambda*I-(kp)*Q*R + km*C) - kappa*Q');
model.component('comp1').physics('tds2').feature('reac2').set('R_Q', '- kappa*Q');
model.component('comp1').physics('tds').feature('cdm1').set('D_I', [0; 0; 0; 0; 0; 0; 0; 0; 0]);
model.component('comp1').physics('tds').feature('init1').set('initc', 'I');
model.component('comp1').physics('tds').feature('reac1').set('R_I', 'mu*C-alpha*I');
model.component('comp1').physics('tds3').feature('cdm1').set('D_R', [0; 0; 0; 0; 0; 0; 0; 0; 0]);
model.component('comp1').physics('tds3').feature('init1').set('initc', 'R');
model.component('comp1').physics('tds3').feature('reac1').set('R_R', 'r - (kp)*Q*R + km*C - beta*R');
model.component('comp1').physics('tds4').feature('cdm1').set('D_C', [0; 0; 0; 0; 0; 0; 0; 0; 0]);
model.component('comp1').physics('tds4').feature('init1').set('initc', 'C');
model.component('comp1').physics('tds4').feature('reac1').set('R_C', 'kp*Q*R - km*C - gamma*C');
if H_3<10
model.component('comp1').mesh('mesh1').feature('size').set('hauto', 2);
model.component('comp1').mesh('mesh1').feature('size').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('ref1').set('boxcoord', true);
model.component('comp1').mesh('mesh1').feature('ref1').set('xmax', 0.1);
model.component('comp1').mesh('mesh1').feature('ref1').set('ymax', 1.1);
model.component('comp1').mesh('mesh1').feature('ref1').set('ymin', 0.9);
model.component('comp1').mesh('mesh1').feature('ref1').set('numrefine', 3);
model.component('comp1').mesh('mesh1').feature('ref2').set('boxcoord', true);
model.component('comp1').mesh('mesh1').feature('ref2').set('xmax', 0.5);
model.component('comp1').mesh('mesh1').feature('ref2').set('ymax', 1.5);
model.component('comp1').mesh('mesh1').feature('ref2').set('ymin', 0.5);
model.component('comp1').mesh('mesh1').feature('ref2').set('numrefine', 2);
model.component('comp1').mesh('mesh1').feature('ref3').set('boxcoord', true);
model.component('comp1').mesh('mesh1').feature('ref3').set('xmin', -0.5);
model.component('comp1').mesh('mesh1').feature('ref3').set('ymax', 1.5);
model.component('comp1').mesh('mesh1').feature('ref3').set('ymin', 1);
if H_3>1
    model.component('comp1').mesh('mesh1').feature('ref4').set('boxcoord', true);
    model.component('comp1').mesh('mesh1').feature('ref4').set('xmin', -5);
    model.component('comp1').mesh('mesh1').feature('ref4').set('xmax', '10+L');
    model.component('comp1').mesh('mesh1').feature('ref4').set('ymax', 2);
    model.component('comp1').mesh('mesh1').feature('ref4').set('numrefine', 2);
end
else
     model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('ref1').set('boxcoord', true);
model.component('comp1').mesh('mesh1').feature('ref1').set('xmax', '5+L');
model.component('comp1').mesh('mesh1').feature('ref1').set('xmin', -5);
model.component('comp1').mesh('mesh1').feature('ref1').set('ymax', 2);
model.component('comp1').mesh('mesh1').feature('ref1').set('numrefine', 1);
model.component('comp1').mesh('mesh1').feature('ref2').set('boxcoord', true);
model.component('comp1').mesh('mesh1').feature('ref2').set('xmax', 'L');
model.component('comp1').mesh('mesh1').feature('ref2').set('ymax', 1.1);
model.component('comp1').mesh('mesh1').feature('ref2').set('ymin', 0.9);
model.component('comp1').mesh('mesh1').feature('ref2').set('numrefine', 2);
end
model.component('comp1').mesh('mesh1').run;

model.study.create('std2');
model.study('std2').create('stat', 'Stationary');
model.study.create('std3');
model.study('std3').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std3');
model.sol('sol1').attach('std3');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol.create('sol2');
model.sol('sol2').study('std2');
model.sol('sol2').attach('std2');
model.sol('sol2').create('st1', 'StudyStep');
model.sol('sol2').create('v1', 'Variables');
model.sol('sol2').create('s1', 'Stationary');
model.sol('sol2').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol2').feature('s1').create('d1', 'Direct');
model.sol('sol2').feature('s1').feature.remove('fcDef');

model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result.create('pg3', 'PlotGroup2D');
model.result.create('pg4', 'PlotGroup2D');
model.result.create('pg5', 'PlotGroup2D');
model.result.create('pg6', 'PlotGroup2D');
model.result.create('pg7', 'PlotGroup2D');
model.result.create('pg8', 'PlotGroup2D');
model.result.create('pg9', 'PlotGroup1D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').create('str1', 'Streamline');
model.result('pg2').create('surf1', 'Surface');
model.result('pg2').create('str1', 'Streamline');
model.result('pg2').feature('surf1').set('expr', 'I');
model.result('pg3').create('surf1', 'Surface');
model.result('pg3').create('str1', 'Streamline');
model.result('pg3').feature('surf1').set('expr', 'R');
model.result('pg4').create('surf1', 'Surface');
model.result('pg4').create('str1', 'Streamline');
model.result('pg4').feature('surf1').set('expr', 'C');
model.result('pg5').set('data', 'dset2');
model.result('pg5').create('surf1', 'Surface');
model.result('pg5').create('str1', 'Streamline');
model.result('pg6').set('data', 'dset2');
model.result('pg6').create('surf1', 'Surface');
model.result('pg6').create('str1', 'Streamline');
model.result('pg6').feature('surf1').set('expr', 'I');
model.result('pg7').set('data', 'dset2');
model.result('pg7').create('surf1', 'Surface');
model.result('pg7').create('str1', 'Streamline');
model.result('pg7').feature('surf1').set('expr', 'R');
model.result('pg8').set('data', 'dset2');
model.result('pg8').create('surf1', 'Surface');
model.result('pg8').create('str1', 'Streamline');
model.result('pg8').feature('surf1').set('expr', 'C');
model.result('pg9').set('data', 'dset2');
model.result('pg9').create('lngr1', 'LineGraph');
model.result('pg9').feature('lngr1').set('xdata', 'expr');
model.result('pg9').feature('lngr1').selection.set([13]);
model.result('pg9').feature('lngr1').set('expr', '-Qy/Q');

model.study('std2').label('Steady QS solver 1');
model.study('std2').feature('stat').set('useinitsol', true);
model.study('std2').feature('stat').set('initmethod', 'sol');
model.study('std2').feature('stat').set('initstudy', 'std3');
model.study('std2').feature('stat').set('solnum', 'last');
model.study('std2').feature('stat').set('usesol', true);
model.study('std2').feature('stat').set('notsolmethod', 'sol');
model.study('std3').label('Unsteady QS solver');
model.study('std3').feature('time').set('tlist', ['range(0,',num2str(t_max),',',num2str(t_max),')']);
model.study('std3').feature('time').set('useinitsol', true);
model.study('std3').feature('time').set('initmethod', 'sol');
model.study('std3').feature('time').set('usesol', true);
model.study('std3').feature('time').set('notsolmethod', 'sol');

model.sol('sol1').attach('std3');
model.sol('sol1').feature('v1').set('initmethod', 'sol');
model.sol('sol1').feature('v1').set('notsolmethod', 'sol');
model.sol('sol1').feature('v1').set('clist', {['range(0,',num2str(t_max),',',num2str(t_max),')'] '10.0[s]'});
model.sol('sol1').feature('t1').set('tlist', ['range(0,',num2str(t_max),',',num2str(t_max),')']);
model.sol('sol1').feature('t1').set('rtol', 0.005);
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 8);
model.sol('sol1').feature('t1').feature('fc1').set('damp', 0.9);
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('t1').feature('fc1').set('stabacc', 'aacc');
model.sol('sol1').feature('t1').feature('fc1').set('aaccdim', 5);
model.sol('sol1').feature('t1').feature('fc1').set('aaccmix', 0.9);
model.sol('sol1').feature('t1').feature('fc1').set('aaccdelay', 1);
model.sol('sol1').feature('t1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').runAll;
model.sol('sol2').attach('std2');
model.sol('sol2').feature('v1').set('initmethod', 'sol');
model.sol('sol2').feature('v1').set('initsol', 'sol1');
model.sol('sol2').feature('v1').set('solnum', 'last');
model.sol('sol2').feature('v1').set('notsolmethod', 'sol');
model.sol('sol2').feature('s1').feature('fc1').set('initstep', 0.01);
model.sol('sol2').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol2').feature('s1').feature('fc1').set('maxiter', 50);
model.sol('sol2').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol2').feature('s1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol2').runAll;

model.result('pg1').label('Concentration (tds2)');
model.result('pg1').set('titletype', 'custom');
model.result('pg1').feature('surf1').set('descr', 'Concentration');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg1').feature('str1').set('posmethod', 'uniform');
model.result('pg1').feature('str1').set('pointtype', 'arrow');
model.result('pg1').feature('str1').set('arrowcount', 30);
model.result('pg1').feature('str1').set('arrowlength', 'logarithmic');
model.result('pg1').feature('str1').set('arrowscale', 78.9744730431727);
model.result('pg1').feature('str1').set('color', 'gray');
model.result('pg1').feature('str1').set('recover', 'pprint');
model.result('pg1').feature('str1').set('arrowcountactive', false);
model.result('pg1').feature('str1').set('arrowscaleactive', false);
model.result('pg1').feature('str1').set('resolution', 'normal');
model.result('pg2').label('Concentration (tds)');
model.result('pg2').set('titletype', 'custom');
model.result('pg2').feature('surf1').set('resolution', 'normal');
model.result('pg2').feature('str1').set('expr', {'tds.tflux_Ix' 'tds.tflux_Iy'});
model.result('pg2').feature('str1').set('posmethod', 'uniform');
model.result('pg2').feature('str1').set('pointtype', 'arrow');
model.result('pg2').feature('str1').set('arrowlength', 'logarithmic');
model.result('pg2').feature('str1').set('color', 'gray');
model.result('pg2').feature('str1').set('recover', 'pprint');
model.result('pg2').feature('str1').set('resolution', 'normal');
model.result('pg3').label('Concentration (tds3)');
model.result('pg3').set('titletype', 'custom');
model.result('pg3').feature('surf1').set('resolution', 'normal');
model.result('pg3').feature('str1').set('expr', {'tds3.tflux_Rx' 'tds3.tflux_Ry'});
model.result('pg3').feature('str1').set('posmethod', 'uniform');
model.result('pg3').feature('str1').set('pointtype', 'arrow');
model.result('pg3').feature('str1').set('arrowlength', 'logarithmic');
model.result('pg3').feature('str1').set('color', 'gray');
model.result('pg3').feature('str1').set('recover', 'pprint');
model.result('pg3').feature('str1').set('resolution', 'normal');
model.result('pg4').label('Concentration (tds4)');
model.result('pg4').set('titletype', 'custom');
model.result('pg4').feature('surf1').set('resolution', 'normal');
model.result('pg4').feature('str1').set('expr', {'tds4.tflux_Cx' 'tds4.tflux_Cy'});
model.result('pg4').feature('str1').set('posmethod', 'uniform');
model.result('pg4').feature('str1').set('pointtype', 'arrow');
model.result('pg4').feature('str1').set('arrowlength', 'logarithmic');
model.result('pg4').feature('str1').set('color', 'gray');
model.result('pg4').feature('str1').set('recover', 'pprint');
model.result('pg4').feature('str1').set('resolution', 'normal');
model.result('pg5').label('Concentration (tds2) 1');
model.result('pg5').set('titletype', 'custom');
model.result('pg5').feature('surf1').set('descr', 'Concentration');
model.result('pg5').feature('surf1').set('resolution', 'normal');
model.result('pg5').feature('str1').set('posmethod', 'uniform');
model.result('pg5').feature('str1').set('pointtype', 'arrow');
model.result('pg5').feature('str1').set('arrowcount', 30);
model.result('pg5').feature('str1').set('arrowlength', 'logarithmic');
model.result('pg5').feature('str1').set('arrowscale', 78.97447304540627);
model.result('pg5').feature('str1').set('color', 'gray');
model.result('pg5').feature('str1').set('recover', 'pprint');
model.result('pg5').feature('str1').set('arrowcountactive', false);
model.result('pg5').feature('str1').set('arrowscaleactive', false);
model.result('pg5').feature('str1').set('resolution', 'normal');
model.result('pg6').label('Concentration (tds) 1');
model.result('pg6').set('titletype', 'custom');
model.result('pg6').feature('surf1').set('resolution', 'normal');
model.result('pg6').feature('str1').set('expr', {'tds.tflux_Ix' 'tds.tflux_Iy'});
model.result('pg6').feature('str1').set('posmethod', 'uniform');
model.result('pg6').feature('str1').set('pointtype', 'arrow');
model.result('pg6').feature('str1').set('arrowlength', 'logarithmic');
model.result('pg6').feature('str1').set('color', 'gray');
model.result('pg6').feature('str1').set('recover', 'pprint');
model.result('pg6').feature('str1').set('resolution', 'normal');
model.result('pg7').label('Concentration (tds3) 1');
model.result('pg7').set('titletype', 'custom');
model.result('pg7').feature('surf1').set('resolution', 'normal');
model.result('pg7').feature('str1').set('expr', {'tds3.tflux_Rx' 'tds3.tflux_Ry'});
model.result('pg7').feature('str1').set('posmethod', 'uniform');
model.result('pg7').feature('str1').set('pointtype', 'arrow');
model.result('pg7').feature('str1').set('arrowlength', 'logarithmic');
model.result('pg7').feature('str1').set('color', 'gray');
model.result('pg7').feature('str1').set('recover', 'pprint');
model.result('pg7').feature('str1').set('resolution', 'normal');
model.result('pg8').label('Concentration (tds4) 1');
model.result('pg8').set('titletype', 'custom');
model.result('pg8').feature('surf1').set('resolution', 'normal');
model.result('pg8').feature('str1').set('expr', {'tds4.tflux_Cx' 'tds4.tflux_Cy'});
model.result('pg8').feature('str1').set('posmethod', 'uniform');
model.result('pg8').feature('str1').set('pointtype', 'arrow');
model.result('pg8').feature('str1').set('arrowlength', 'logarithmic');
model.result('pg8').feature('str1').set('color', 'gray');
model.result('pg8').feature('str1').set('recover', 'pprint');
model.result('pg8').feature('str1').set('resolution', 'normal');

model.result.numerical.create('max1', 'MaxSurface');
model.result.numerical('max1').selection.set([2]);
model.result.numerical('max1').setIndex('expr', 'Q', 0);
model.result.numerical('max1').set('data', 'dset2');
model.result.table.create('tbl1', 'Table');
model.result.table('tbl1').comments('Surface Maximum 1');
model.result.numerical('max1').set('table', 'tbl1');
model.result.numerical('max1').setResult;

model.result.numerical.create('min1', 'MinSurface');
model.result.numerical('min1').selection.set([2]);
model.result.numerical('min1').setIndex('expr', 'Q', 0);
model.result.numerical('min1').set('data', 'dset2');
model.result.table.create('tbl2', 'Table');
model.result.table('tbl2').comments('Surface Minimum 1');
model.result.numerical('min1').set('table', 'tbl2');
model.result.numerical('min1').setResult;

model.result.numerical.create('max2', 'MaxSurface');
model.result.numerical('max2').selection.set([2]);
model.result.numerical('max2').setIndex('expr', 'C', 0);
model.result.numerical('max2').set('data', 'dset2');
model.result.table.create('tbl3', 'Table');
model.result.table('tbl3').comments('Surface Maximum 2');
model.result.numerical('max2').set('table', 'tbl3');
model.result.numerical('max2').setResult;

model.result.numerical.create('pev1', 'EvalPoint');
model.result.numerical('pev1').set('data', 'dset2');
model.result.numerical('pev1').selection.set([3]);
model.result.numerical('pev1').set('probetag', 'none');
model.result.table.create('tbl4', 'Table');
model.result.table('tbl4').comments('Bottom Left Q');
model.result.numerical('pev1').set('table', 'tbl4');
model.result.numerical('pev1').set('expr', {'Q'});
model.result.numerical('pev1').set('unit', {''});
model.result.numerical('pev1').set('descr', {'Concentration'});
model.result.numerical('pev1').setResult;


model.result('pg9').set('xlabel', 'x-coordinate');
model.result('pg9').set('ylabel', '-Qy/Q');
model.result('pg9').set('xlabelactive', false);
model.result('pg9').set('ylabelactive', false);
model.result('pg9').feature('lngr1').set('xdataexpr', 'x');
model.result('pg9').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg9').feature('lngr1').set('resolution', 'coarse');
model.result('pg9').feature('lngr1').set('smooth', 'everywhere');
model.result('pg9').feature('lngr1').set('recover', 'ppr');
model.result('pg9').feature('lngr1').set('resolution', 'coarse');
figure
pd = mphplot(model,'pg9');
close
x = pd{1}{1}.p;
Pe_eff = pd{1}{1}.d';

out = model;
Qmax = model.result.table('tbl1').getReal;
Qmin = model.result.table('tbl2').getReal;
Cmax = model.result.table('tbl3').getReal;
Qleft = model.result.table('tbl4').getReal;