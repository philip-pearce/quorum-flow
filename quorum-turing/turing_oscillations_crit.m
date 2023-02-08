function [Amin,Amax,uxav,uyav,Tc] = turing_oscillations_crit(Tsweep,Ksweep,kfa)
%
% turing_oscillations.m
%
% Model exported on Nov 10 2021, 10:15 by COMSOL 5.6.0.401.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');
model.param.set('fu', '0.28');
model.param.set('fv', '-0.5');
model.param.set('gu', '0.75');
model.param.set('gv', '-0.5');
model.param.set('du', '70');
model.param.set('dv', '875');
model.param.set('L', '717');

model.component.create('comp3', true);

model.component('comp3').geom.create('geom3', 2);

model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');
model.result.table.create('tbl4', 'Table');
model.result.table.create('tbl5', 'Table');

model.component('comp3').mesh.create('mesh3');

model.component('comp3').geom('geom3').create('sq1', 'Square');
model.component('comp3').geom('geom3').feature('sq1').set('size', 717);
model.component('comp3').geom('geom3').run;

model.frame('material1').tag('material3');
model.frame('mesh1').tag('mesh3');
model.frame('geometry1').tag('geometry3');
model.frame('spatial1').tag('spatial3');

model.variable.create('var1');

model.component('comp3').view('view1').tag('view4');
model.view.create('view3', 2);

model.component('comp3').physics.create('tds2', 'DilutedSpecies', 'geom3');
model.component('comp3').physics('tds2').identifier('tds2');
model.component('comp3').physics('tds2').field('concentration').field('u');
model.component('comp3').physics('tds2').field('concentration').component({'u'});
model.component('comp3').physics('tds2').create('reac1', 'Reactions', 2);
model.component('comp3').physics('tds2').feature('reac1').selection.set([1]);
model.component('comp3').physics.create('tds5', 'DilutedSpecies', 'geom3');
model.component('comp3').physics('tds5').identifier('tds5');
model.component('comp3').physics('tds5').field('concentration').field('v');
model.component('comp3').physics('tds5').field('concentration').component({'v'});
model.component('comp3').physics('tds5').create('reac1', 'Reactions', 2);
model.component('comp3').physics('tds5').feature('reac1').selection.set([1]);

model.component('comp3').mesh('mesh3').create('ftri1', 'FreeTri');

model.result.table('tbl1').comments('Global Evaluation 1');
model.result.table('tbl2').comments('Line Maximum 1');
model.result.table('tbl3').comments('Surface Maximum 1');
model.result.table('tbl4').comments('Surface Average 1');
model.result.table('tbl5').comments('Surface Average 2');

model.view('view3').axis.set('xmin', -0.04999998211860657);
model.view('view3').axis.set('xmax', 1.0499999523162842);
model.view('view3').axis.set('ymin', -0.773396372795105);
model.view('view3').axis.set('ymax', 0.8733963966369629);
model.component('comp3').view('view4').label('View 4');
model.component('comp3').view('view4').axis.set('xmin', -123.80267333984375);
model.component('comp3').view('view4').axis.set('xmax', 840.8026733398438);
model.component('comp3').view('view4').axis.set('ymin', -147.37173461914062);
model.component('comp3').view('view4').axis.set('ymax', 864.3717041015625);

model.component('comp3').physics('tds2').label('Transport of Diluted Species');
model.component('comp3').physics('tds2').feature('cdm1').set('D_u', {'du'; '0'; '0'; '0'; 'du'; '0'; '0'; '0'; 'du'});
model.component('comp3').physics('tds2').feature('init1').set('initc', '0.001*sin(pi*y/100)');
model.component('comp3').physics('tds2').feature('reac1').set('R_u', '(fu + kfu * x/717) * u + fv * v - u^3');
model.component('comp3').physics('tds5').label('Transport of Diluted Species 2');
model.component('comp3').physics('tds5').feature('cdm1').set('D_v', {'dv'; '0'; '0'; '0'; 'dv'; '0'; '0'; '0'; 'dv'});
model.component('comp3').physics('tds5').feature('reac1').set('R_v', 'gu * u +gv * v');

model.component('comp3').mesh('mesh3').feature('size').set('hauto', 1);
model.component('comp3').mesh('mesh3').feature('size').set('table', 'cfd');
model.component('comp3').mesh('mesh3').run;

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');
model.study.create('std2');
model.study('std2').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol.create('sol2');
model.sol('sol2').study('std2');
model.sol('sol2').attach('std2');
model.sol('sol2').create('st1', 'StudyStep');
model.sol('sol2').create('v1', 'Variables');
model.sol('sol2').create('t1', 'Time');
model.sol('sol2').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol2').feature('t1').create('d1', 'Direct');
model.sol('sol2').feature('t1').create('i1', 'Iterative');
model.sol('sol2').feature('t1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol2').feature('t1').feature.remove('fcDef');

model.result.dataset.create('cln1', 'CutLine2D');
model.result.dataset.create('dset3', 'Solution');
model.result.dataset.create('dset4', 'Solution');
model.result.dataset.create('dset5', 'Solution');
model.result.dataset.create('dset6', 'Solution');
model.result.dataset.create('dset7', 'Solution');
model.result.dataset('cln1').set('data', 'none');
model.result.dataset('dset3').set('solution', 'none');
model.result.dataset('dset4').set('solution', 'none');
model.result.dataset('dset5').set('solution', 'none');
model.result.dataset('dset6').set('solution', 'none');
model.result.dataset('dset7').set('solution', 'none');
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical.create('max1', 'MaxSurface');
model.result.numerical.create('min1', 'MinSurface');
model.result.numerical.create('av1', 'AvSurface');
model.result.numerical.create('av2', 'AvSurface');
model.result.numerical('gev1').set('data', 'dset2');
model.result.numerical('gev1').set('probetag', 'none');
model.result.numerical('max1').set('data', 'dset2');
model.result.numerical('max1').selection.set([1]);
model.result.numerical('max1').set('probetag', 'none');
model.result.numerical('min1').set('data', 'dset2');
model.result.numerical('min1').selection.set([1]);
model.result.numerical('min1').set('probetag', 'none');

model.result.numerical('av1').set('data', 'dset2');
model.result.numerical('av1').selection.set([1]);
model.result.numerical('av1').set('probetag', 'none');

model.result.numerical('av2').set('data', 'dset2');
model.result.numerical('av2').selection.set([1]);
model.result.numerical('av2').set('probetag', 'none');

model.result.create('pg5', 'PlotGroup1D');
model.result.create('pg6', 'PlotGroup1D');
model.result.create('pg7', 'PlotGroup1D');
model.result.create('pg8', 'PlotGroup1D');
model.result.create('pg9', 'PlotGroup1D');
model.result.create('pg10', 'PlotGroup1D');
model.result.create('pg11', 'PlotGroup1D');
model.result.create('pg12', 'PlotGroup1D');
model.result.create('pg13', 'PlotGroup1D');
model.result.create('pg14', 'PlotGroup1D');
model.result.create('pg15', 'PlotGroup1D');
model.result.create('pg16', 'PlotGroup2D');
model.result('pg5').set('data', 'dset2');
model.result('pg5').create('lngr1', 'LineGraph');
model.result('pg5').create('lngr2', 'LineGraph');
model.result('pg5').feature('lngr1').set('xdata', 'expr');
model.result('pg5').feature('lngr2').set('xdata', 'expr');
model.result('pg5').feature('lngr2').set('expr', 'beta>betamin');
model.result('pg6').create('lngr1', 'LineGraph');
model.result('pg6').create('lngr2', 'LineGraph');
model.result('pg6').feature('lngr1').set('xdata', 'expr');
model.result('pg6').feature('lngr2').set('expr', '0.1*((fu + kfu * x/717)>0.3)');
model.result('pg7').create('lngr1', 'LineGraph');
model.result('pg7').feature('lngr1').set('xdata', 'expr');
model.result('pg7').feature('lngr1').set('expr', 'tds4.v');
model.result('pg8').create('tblp1', 'Table');
model.result('pg9').create('tblp1', 'Table');
model.result('pg10').set('data', 'dset5');
model.result('pg10').create('lngr1', 'LineGraph');
model.result('pg10').feature('lngr1').set('xdata', 'expr');
model.result('pg10').feature('lngr1').set('expr', 'tds3.u');
model.result('pg11').set('data', 'dset5');
model.result('pg11').create('lngr1', 'LineGraph');
model.result('pg11').feature('lngr1').set('xdata', 'expr');
model.result('pg11').feature('lngr1').set('expr', 'tds4.v');
model.result('pg12').set('data', 'dset6');
model.result('pg12').create('lngr1', 'LineGraph');
model.result('pg12').feature('lngr1').set('xdata', 'expr');
model.result('pg12').feature('lngr1').set('expr', 'tds3.u');
model.result('pg13').set('data', 'dset6');
model.result('pg13').create('lngr1', 'LineGraph');
model.result('pg13').feature('lngr1').set('xdata', 'expr');
model.result('pg13').feature('lngr1').set('expr', 'tds4.v');
model.result('pg14').set('data', 'dset7');
model.result('pg14').create('lngr1', 'LineGraph');
model.result('pg14').feature('lngr1').set('xdata', 'expr');
model.result('pg14').feature('lngr1').set('expr', 'tds3.u');
model.result('pg15').set('data', 'dset7');
model.result('pg15').create('lngr1', 'LineGraph');
model.result('pg15').feature('lngr1').set('xdata', 'expr');
model.result('pg15').feature('lngr1').set('expr', 'tds4.v');
model.result('pg16').set('data', 'dset2');
model.result('pg16').create('surf1', 'Surface');

model.study('std1').feature('stat').set('useinitsol', true);
model.study('std1').feature('stat').set('initmethod', 'sol');
model.study('std1').feature('stat').set('initstudy', 'std2');
model.study('std1').feature('stat').set('solnum', 'last');
model.study('std1').feature('stat').set('auto_ngen', 15);
model.study('std1').feature('stat').set('manual_ngen', 15);
model.study('std1').feature('stat').set('auto_ngenactive', false);
model.study('std1').feature('stat').set('manual_ngenactive', false);
model.study('std2').feature('time').set('usertol', true);
model.study('std2').feature('time').set('rtol', '0.00005');
model.study('std2').feature('time').set('plot', true);
model.study('std2').feature('time').set('plotgroup', 'pg16');
model.study('std2').feature('time').set('plotfreq', 'tsteps');
model.study('std2').feature('time').set('useinitsol', true);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Stationary');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('initmethod', 'sol');
model.sol('sol1').feature('v1').set('initsol', 'sol2');
model.sol('sol1').feature('v1').set('solnum', 'last');
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'd1');
model.sol('sol1').feature('s1').feature('fc1').set('dtech', 'hnlin');
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 50);
model.sol('sol1').feature('s1').feature('d1').label('Direct, concentrations (tds3) (merged)');
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s1').feature('i1').label('AMG, concentrations (tds4)');
model.sol('sol1').feature('s1').feature('i1').set('nlinnormuse', true);
model.sol('sol1').feature('s1').feature('i1').set('maxlinit', 1000);
model.sol('sol1').feature('s1').feature('i1').feature('ilDef').label('Incomplete LU 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').label('Multigrid 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('maxcoarsedof', 50000);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').label('SOR Line 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('iter', 1);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.7);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('relax', 0.5);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').label('SOR Line 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('iter', 1);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.7);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('relax', 0.5);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);

model.sol('sol2').attach('std2');
model.sol('sol2').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol2').feature('v1').label('Dependent Variables 1.1');

model.sol('sol2').feature('v1').set('scalemethod', 'manual');
model.sol('sol2').feature('v1').set('scaleval', 0.1);
model.sol('sol2').feature('t1').label('Time-Dependent Solver 1.1');

model.sol('sol2').feature('t1').set('rtol', '0.00005');
model.sol('sol2').feature('t1').set('atolglobalfactor', 0.01);

model.sol('sol2').feature('t1').set('maxorder', 2);
model.sol('sol2').feature('t1').set('stabcntrl', true);
model.sol('sol2').feature('t1').set('plot', true);
model.sol('sol2').feature('t1').set('plotgroup', 'pg16');
model.sol('sol2').feature('t1').set('plotfreq', 'tsteps');
model.sol('sol2').feature('t1').feature('dDef').label('Direct 2');
model.sol('sol2').feature('t1').feature('aDef').label('Advanced 1');
model.sol('sol2').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol2').feature('t1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol2').feature('t1').feature('fc1').set('linsolver', 'd1');
model.sol('sol2').feature('t1').feature('fc1').set('maxiter', 8);
model.sol('sol2').feature('t1').feature('fc1').set('damp', 0.9);
model.sol('sol2').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol2').feature('t1').feature('fc1').set('stabacc', 'aacc');
model.sol('sol2').feature('t1').feature('fc1').set('aaccdim', 5);
model.sol('sol2').feature('t1').feature('fc1').set('aaccmix', 0.9);
model.sol('sol2').feature('t1').feature('fc1').set('aaccdelay', 1);
model.sol('sol2').feature('t1').feature('d1').label('Direct, concentrations (tds3) (merged)');
model.sol('sol2').feature('t1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol2').feature('t1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol2').feature('t1').feature('i1').label('AMG, concentrations (tds4)');
model.sol('sol2').feature('t1').feature('i1').set('maxlinit', 50);
model.sol('sol2').feature('t1').feature('i1').feature('ilDef').label('Incomplete LU 1');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').label('Multigrid 1.1');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').set('maxcoarsedof', 50000);
model.sol('sol2').feature('t1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol2').feature('t1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').label('SOR Line 1.1');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('iter', 1);
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.7);
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('relax', 0.5);
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').label('SOR Line 1.1');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('iter', 1);
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.7);
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('relax', 0.5);
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol2').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);


model.variable('var1').set('kfu', num2str(kfa));
model.study('std2').feature('time').set('initmethod', 'init');
model.study('std2').feature('time').set('initstudy', 'zero');
model.sol('sol2').feature('v1').set('clist', {'range(0,1e4,1e4)' '100.0[s]'});
model.sol('sol2').feature('t1').set('tlist', 'range(0,1e4,1e4)');
model.sol('sol2').feature('t1').set('maxstepconstraintbdf', 'auto');
disp('Performing initial time-dependent solve...')
model.sol('sol2').runAll;

disp('Performing initial steady solve...')
model.sol('sol1').runAll;

for T = Tsweep
    for K=Ksweep
model.param.set('T', num2str(T));
model.param.set('k', num2str(K));
model.variable('var1').set('kfu', [num2str(kfa),'*(1 + k*sin(2*pi*t/T))']);
model.sol('sol2').feature('v1').set('clist', {'range(0,T/100,T*10)' '100.0[s]'});
model.sol('sol2').feature('t1').set('tlist', 'range(0,T/100,T*10)');
model.sol('sol2').feature('t1').set('maxstepconstraintbdf', 'const');
model.sol('sol2').feature('t1').set('maxstepbdf', 'T/10');
model.sol('sol2').feature('v1').set('initmethod', 'sol');
model.sol('sol2').feature('v1').set('initsol', 'sol1');
model.sol('sol2').feature('v1').set('solnum', 'last');
disp('Performing full time-dependent solve...')
model.study('std2').feature('time').set('rtol', '0.00005');
model.sol('sol2').runAll;

model.result.dataset('cln1').set('genpoints', [0 0.05; 1 0.05]);
model.result.numerical('gev1').set('table', 'tbl1');
model.result.numerical('gev1').set('expr', {'kfu'});
model.result.numerical('gev1').set('unit', {''});
model.result.numerical('gev1').set('descr', {''});
model.result.numerical('max1').set('expr', {'u'});
model.result.numerical('max1').set('unit', {'mol/m^3'});
model.result.numerical('max1').set('descr', {'Concentration'});
model.result.numerical('min1').set('expr', {'u'});
model.result.numerical('min1').set('unit', {'mol/m^3'});
model.result.numerical('min1').set('descr', {'Concentration'});

model.result.numerical('av1').set('expr', {'abs(ux)'});
model.result.numerical('av1').set('unit', {'mol/m^4'});
model.result.numerical('av1').set('descr', {'Concentration'});

model.result.numerical('av1').set('expr', {'abs(uy)'});
model.result.numerical('av1').set('unit', {'mol/m^4'});
model.result.numerical('av1').set('descr', {'Concentration'});

model.result.numerical('gev1').setResult;
model.result('pg5').set('looplevelinput', {'manual'});
model.result('pg5').set('looplevel', [1]);
model.result('pg5').set('xlabel', 'Spatial coordinates, x component (m)');
model.result('pg5').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg5').set('xlabelactive', false);
model.result('pg5').set('ylabelactive', false);
model.result('pg5').feature('lngr1').set('descr', 'Concentration');
model.result('pg5').feature('lngr1').set('xdataexpr', 'x');
model.result('pg5').feature('lngr1').set('xdataunit', 'm');
model.result('pg5').feature('lngr1').set('xdatadescr', 'Spatial coordinates, x component');
model.result('pg5').feature('lngr1').set('smooth', 'internal');
model.result('pg5').feature('lngr1').set('resolution', 'normal');
model.result('pg5').feature('lngr2').active(false);
model.result('pg5').feature('lngr2').set('xdataexpr', 'x');
model.result('pg5').feature('lngr2').set('xdataunit', 'm');
model.result('pg5').feature('lngr2').set('xdatadescr', 'Spatial coordinates, x component');
model.result('pg5').feature('lngr2').set('smooth', 'internal');
model.result('pg5').feature('lngr2').set('resolution', 'normal');
model.result('pg6').label('Concentration (tds3)');
model.result('pg6').set('titletype', 'custom');
model.result('pg6').set('typeintitle', false);
model.result('pg6').set('xlabel', 'Spatial coordinates, x component (m)');
model.result('pg6').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg6').set('xlabelactive', false);
model.result('pg6').set('ylabelactive', false);
model.result('pg6').feature('lngr1').set('descr', 'Concentration');
model.result('pg6').feature('lngr1').set('xdataexpr', 'x');
model.result('pg6').feature('lngr1').set('xdataunit', 'm');
model.result('pg6').feature('lngr1').set('xdatadescr', 'Spatial coordinates, x component');
model.result('pg6').feature('lngr1').set('smooth', 'internal');
model.result('pg6').feature('lngr1').set('resolution', 'normal');
model.result('pg6').feature('lngr2').active(false);
model.result('pg6').feature('lngr2').set('smooth', 'internal');
model.result('pg6').feature('lngr2').set('resolution', 'normal');
model.result('pg7').label('Concentration (tds4)');
model.result('pg7').set('titletype', 'custom');
model.result('pg7').set('typeintitle', false);
model.result('pg7').set('xlabel', 'Spatial coordinates, x component (m)');
model.result('pg7').set('ylabel', 'Velocity field, y component (m/s)');
model.result('pg7').set('xlabelactive', false);
model.result('pg7').set('ylabelactive', false);
model.result('pg7').feature('lngr1').set('descr', 'Velocity field, y component');
model.result('pg7').feature('lngr1').set('xdataexpr', 'x');
model.result('pg7').feature('lngr1').set('xdataunit', 'm');
model.result('pg7').feature('lngr1').set('xdatadescr', 'Spatial coordinates, x component');
model.result('pg7').feature('lngr1').set('smooth', 'internal');
model.result('pg7').feature('lngr1').set('resolution', 'normal');
model.result('pg8').set('data', 'none');
model.result('pg8').set('xlabel', 'Time (s)');
model.result('pg8').set('ylabel', 'kfu');
model.result('pg8').set('xlabelactive', false);
model.result('pg8').set('ylabelactive', false);
model.result('pg9').set('data', 'none');
model.result('pg9').set('xlabel', 'Time (s)');
model.result('pg9').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg9').set('xlabelactive', false);
model.result('pg9').set('ylabelactive', false);
model.result('pg9').feature('tblp1').set('table', 'tbl2');
model.result('pg10').label('Concentration (tds3) 1');
model.result('pg10').set('titletype', 'custom');
model.result('pg10').set('typeintitle', false);
model.result('pg10').set('xlabel', 'x-coordinate (m)');
model.result('pg10').set('ylabel', 'Velocity field, x component');
model.result('pg10').set('xlabelactive', false);
model.result('pg10').set('ylabelactive', false);
model.result('pg10').feature('lngr1').set('descr', 'Velocity field, x component');
model.result('pg10').feature('lngr1').set('xdataexpr', 'x');
model.result('pg10').feature('lngr1').set('xdataunit', 'm');
model.result('pg10').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg10').feature('lngr1').set('smooth', 'internal');
model.result('pg10').feature('lngr1').set('allowmaterialsmoothing', false);
model.result('pg10').feature('lngr1').set('resolution', 'normal');
model.result('pg11').label('Concentration (tds4) 1');
model.result('pg11').set('titletype', 'custom');
model.result('pg11').set('typeintitle', false);
model.result('pg11').feature('lngr1').set('unit', 'm/s');
model.result('pg11').feature('lngr1').set('descr', 'Velocity field, y component');
model.result('pg11').feature('lngr1').set('xdataexpr', 'x');
model.result('pg11').feature('lngr1').set('xdataunit', 'm');
model.result('pg11').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg11').feature('lngr1').set('smooth', 'internal');
model.result('pg11').feature('lngr1').set('allowmaterialsmoothing', false);
model.result('pg11').feature('lngr1').set('resolution', 'normal');
model.result('pg12').label('Concentration (tds3) 2');
model.result('pg12').set('titletype', 'custom');
model.result('pg12').set('typeintitle', false);
model.result('pg12').set('xlabel', 'x-coordinate (m)');
model.result('pg12').set('ylabel', 'Velocity field, x component (m/s)');
model.result('pg12').set('xlabelactive', false);
model.result('pg12').set('ylabelactive', false);
model.result('pg12').feature('lngr1').set('unit', 'm/s');
model.result('pg12').feature('lngr1').set('descr', 'Velocity field, x component');
model.result('pg12').feature('lngr1').set('xdataexpr', 'x');
model.result('pg12').feature('lngr1').set('xdataunit', 'm');
model.result('pg12').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg12').feature('lngr1').set('smooth', 'internal');
model.result('pg12').feature('lngr1').set('allowmaterialsmoothing', false);
model.result('pg12').feature('lngr1').set('resolution', 'normal');
model.result('pg13').label('Concentration (tds4) 2');
model.result('pg13').set('titletype', 'custom');
model.result('pg13').set('typeintitle', false);
model.result('pg13').feature('lngr1').set('unit', 'm/s');
model.result('pg13').feature('lngr1').set('descr', 'Velocity field, y component');
model.result('pg13').feature('lngr1').set('xdataexpr', 'x');
model.result('pg13').feature('lngr1').set('xdataunit', 'm');
model.result('pg13').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg13').feature('lngr1').set('smooth', 'internal');
model.result('pg13').feature('lngr1').set('allowmaterialsmoothing', false);
model.result('pg13').feature('lngr1').set('resolution', 'normal');
model.result('pg14').label('Concentration (tds3) 3');
model.result('pg14').set('titletype', 'custom');
model.result('pg14').set('typeintitle', false);
model.result('pg14').set('xlabel', 'x-coordinate (m)');
model.result('pg14').set('ylabel', 'Velocity field, x component (m/s)');
model.result('pg14').set('xlabelactive', false);
model.result('pg14').set('ylabelactive', false);
model.result('pg14').feature('lngr1').set('unit', 'm/s');
model.result('pg14').feature('lngr1').set('descr', 'Velocity field, x component');
model.result('pg14').feature('lngr1').set('xdataexpr', 'x');
model.result('pg14').feature('lngr1').set('xdataunit', 'm');
model.result('pg14').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg14').feature('lngr1').set('smooth', 'internal');
model.result('pg14').feature('lngr1').set('allowmaterialsmoothing', false);
model.result('pg14').feature('lngr1').set('resolution', 'normal');
model.result('pg15').label('Concentration (tds4) 3');
model.result('pg15').set('titletype', 'custom');
model.result('pg15').set('typeintitle', false);
model.result('pg15').feature('lngr1').set('unit', 'm/s');
model.result('pg15').feature('lngr1').set('descr', 'Velocity field, y component');
model.result('pg15').feature('lngr1').set('xdataexpr', 'x');
model.result('pg15').feature('lngr1').set('xdataunit', 'm');
model.result('pg15').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg15').feature('lngr1').set('smooth', 'internal');
model.result('pg15').feature('lngr1').set('allowmaterialsmoothing', false);
model.result('pg15').feature('lngr1').set('resolution', 'normal');
model.result('pg16').set('looplevel', [101]);
model.result('pg16').feature('surf1').set('resolution', 'normal');

model.result.numerical('max1').set('table', 'tbl3');
model.result.numerical('max1').set('data', 'dset2');
model.result.numerical('max1').set('table', 'tbl3');
model.result.numerical('max1').set('expr', {'u'});
model.result.numerical('max1').set('unit', {''});
model.result.numerical('max1').set('descr', {''});
model.result.numerical('max1').setResult;
model.result.numerical('min1').set('table', 'tbl2');
model.result.numerical('min1').set('data', 'dset2');
model.result.numerical('min1').set('table', 'tbl2');
model.result.numerical('min1').set('expr', {'u'});
model.result.numerical('min1').set('unit', {''});
model.result.numerical('min1').set('descr', {''});
model.result.numerical('min1').setResult;

model.result.numerical('av1').set('table', 'tbl4');
model.result.numerical('av1').set('data', 'dset2');
model.result.numerical('av1').set('table', 'tbl4');
model.result.numerical('av1').set('expr', {'abs(ux)'});
model.result.numerical('av1').set('unit', {''});
model.result.numerical('av1').set('descr', {''});
model.result.numerical('av1').setResult;

model.result.numerical('av2').set('table', 'tbl5');
model.result.numerical('av2').set('data', 'dset2');
model.result.numerical('av2').set('table', 'tbl5');
model.result.numerical('av2').set('expr', {'abs(uy)'});
model.result.numerical('av2').set('unit', {''});
model.result.numerical('av2').set('descr', {''});
model.result.numerical('av2').setResult;

Amax = model.result.table('tbl3').getReal;
Amin = model.result.table('tbl2').getReal;
uxav = model.result.table('tbl4').getReal;
uyav = model.result.table('tbl5').getReal;

save(['output_crit_T_',num2str(T),'_K_',num2str(K),'_kfa_',num2str(kfa),'.mat'],'Amax','Amin','uxav','uyav')
flag = 0;

if min(Amax/(0.5^(1/3) * 875^(1/6) / (717^(1/3))))<2*pi*sqrt(875 / (0.5*717^2))/sqrt(sqrt(0.5*0.75/0.5^2/(70/875))-1)
flag = 1;
    break
end
    
    end
if flag==1
break
end
end

Tc = T;
