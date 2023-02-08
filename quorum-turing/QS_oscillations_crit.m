function [Qmin,Qmax,Tc] = QS_oscillations_crit(Tsweep,Ksweep,kA)
%
% QS_oscillations.m
%
% Model exported on Nov 9 2021, 13:23 by COMSOL 5.6.0.401.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.label('QS_oscillations.mph');

model.baseSystem('none');

model.param.set('q', '0.04');
model.param.set('r', '4e-4');
model.param.set('lambda', '0.002');
model.param.set('mu', '200');
model.param.set('kpb', '0');
model.param.set('km', '0.002');
model.param.set('alpha', '0.0002');
model.param.set('beta', '0.002');
model.param.set('gamma', '0.0002');
model.param.set('kappa', '2e-05');
model.param.set('Pe_mean', '2000');
model.param.set('L', '100'); %note set mesh hauto=2 if L=2, hauto=3 if L=50
model.param.set('t_max', '20/epsilon');
model.param.set('epsilon', '2e-6');
model.param.set('rho', '0.1');

model.param.set('kA', num2str(kA));

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');
model.result.table.create('tbl4', 'Table');
model.result.table.create('tbl5', 'Table');
model.result.table.create('tbl6', 'Table');

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').label('cells');
model.component('comp1').geom('geom1').feature('r2').set('pos', [0 0]);
model.component('comp1').geom('geom1').feature('r2').set('size', {'L' '1'});
model.component('comp1').geom('geom1').create('pard1', 'PartitionDomains');
model.component('comp1').geom('geom1').feature('pard1').set('partitionwith', 'edges');
model.component('comp1').geom('geom1').feature('pard1').selection('domain').set('r2(1)', 1);
model.component('comp1').geom('geom1').feature('pard1').selection('edge').set('r2(1)', 1);
model.component('comp1').geom('geom1').run;

model.variable.create('var1');


model.component('comp1').physics.create('tds2', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds2').identifier('tds2');
model.component('comp1').physics('tds2').field('concentration').field('Q');
model.component('comp1').physics('tds2').field('concentration').component({'Q'});
model.component('comp1').physics('tds2').create('reac1', 'Reactions', 2);
model.component('comp1').physics('tds2').feature('reac1').selection.set([1]);
model.component('comp1').physics('tds2').create('fl1', 'FluxBoundary', 1);
model.component('comp1').physics('tds2').feature('fl1').selection.set([3]);
model.component('comp1').physics.create('tds', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds').field('concentration').field('I');
model.component('comp1').physics('tds').field('concentration').component({'I'});
model.component('comp1').physics('tds').create('reac1', 'Reactions', 2);
model.component('comp1').physics('tds').feature('reac1').selection.set([1]);
model.component('comp1').physics.create('tds3', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds3').field('concentration').field('R');
model.component('comp1').physics('tds3').field('concentration').component({'R'});
model.component('comp1').physics('tds3').create('reac1', 'Reactions', 2);
model.component('comp1').physics('tds3').feature('reac1').selection.set([1]);
model.component('comp1').physics.create('tds4', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds4').field('concentration').field('C');
model.component('comp1').physics('tds4').field('concentration').component({'C'});
model.component('comp1').physics('tds4').create('reac1', 'Reactions', 2);
model.component('comp1').physics('tds4').feature('reac1').selection.set([1]);

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');

model.component('comp1').probe.create('pdom1', 'DomainPoint');

model.result.table('tbl1').comments('Surface Maximum 1');
model.result.table('tbl2').comments('Point Evaluation 1');
model.result.table('tbl3').label('Probe Table 3');
model.result.table('tbl4').comments('Point Probe Expression 1');
model.result.table('tbl5').comments('Surface Minimum 1');
model.result.table('tbl6').comments('Surface Minimum 1');


model.component('comp1').view('view1').axis.set('xmin', -0.04999999701976776);
model.component('comp1').view('view1').axis.set('xmax', 2.049999952316284);
model.component('comp1').view('view1').axis.set('ymin', -0.7933672666549683);
model.component('comp1').view('view1').axis.set('ymax', 1.7933672666549683);

model.component('comp1').physics('tds2').feature('cdm1').set('D_Q', [1; 0; 0; 0; 1; 0; 0; 0; 1]);
model.component('comp1').physics('tds2').feature('init1').set('initc', 'Q');
model.component('comp1').physics('tds2').feature('reac1').set('R_Q', 'rho*(q+lambda*I-(kp)*Q*R + km*C) - kappa*Q');
model.component('comp1').physics('tds2').feature('fl1').set('FluxType', 'ExternalConvection');
model.component('comp1').physics('tds2').feature('fl1').set('species', true);
model.component('comp1').physics('tds2').feature('fl1').set('kc', 1);
model.component('comp1').physics('tds').feature('cdm1').set('D_I', [0; 0; 0; 0; 0; 0; 0; 0; 0]);
model.component('comp1').physics('tds').feature('init1').set('initc', 'I');
model.component('comp1').physics('tds').feature('reac1').set('R_I', 'mu*C-alpha*I');
model.component('comp1').physics('tds3').feature('cdm1').set('D_R', [0; 0; 0; 0; 0; 0; 0; 0; 0]);
model.component('comp1').physics('tds3').feature('init1').set('initc', 'R');
model.component('comp1').physics('tds3').feature('reac1').set('R_R', 'r - (kp)*Q*R + km*C - beta*R');
model.component('comp1').physics('tds4').feature('cdm1').set('D_C', [0; 0; 0; 0; 0; 0; 0; 0; 0]);
model.component('comp1').physics('tds4').feature('init1').set('initc', 'C');
model.component('comp1').physics('tds4').feature('reac1').set('R_C', 'kp*Q*R - km*C - gamma*C');

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 3);
model.component('comp1').mesh('mesh1').feature('size').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('ftri1').set('smoothmaxiter', 8);
model.component('comp1').mesh('mesh1').feature('ftri1').set('smoothmaxdepth', 8);
model.component('comp1').mesh('mesh1').run;

model.component('comp1').probe('pdom1').set('coords2', [2 0]);
model.component('comp1').probe('pdom1').feature('ppb1').set('table', 'tbl3');
model.component('comp1').probe('pdom1').feature('ppb1').set('window', 'window1');

model.study.create('std2');
model.study('std2').create('stat', 'Stationary');
model.study.create('std3');
model.study('std3').create('time', 'Transient');

model.sol.create('sol2');
model.sol('sol2').study('std2');
model.sol('sol2').attach('std2');
model.sol('sol2').create('st1', 'StudyStep');
model.sol('sol2').create('v1', 'Variables');
model.sol('sol2').create('s1', 'Stationary');
model.sol('sol2').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol2').feature('s1').create('d1', 'Direct');
model.sol('sol2').feature('s1').feature.remove('fcDef');
model.sol.create('sol3');
model.sol('sol3').study('std3');
model.sol('sol3').attach('std3');
model.sol('sol3').create('st1', 'StudyStep');
model.sol('sol3').create('v1', 'Variables');
model.sol('sol3').create('t1', 'Time');
model.sol('sol3').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol3').feature('t1').create('d1', 'Direct');
model.sol('sol3').feature('t1').create('i1', 'Iterative');
model.sol('sol3').feature('t1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol3').feature('t1').feature.remove('fcDef');

model.result.dataset.create('dset3', 'Solution');
model.result.dataset.create('cpt1', 'CutPoint2D');
model.result.dataset.create('dset4', 'Solution');
model.result.dataset.create('dset5', 'Solution');
model.result.dataset.create('dset6', 'Solution');
model.result.dataset.create('dset7', 'Solution');
model.result.dataset.create('dset8', 'Solution');
model.result.dataset.create('dset9', 'Solution');
model.result.dataset.create('dset10', 'Solution');
model.result.dataset('dset1').set('solution', 'none');
model.result.dataset('dset2').set('solution', 'sol2');
model.result.dataset('dset3').set('probetag', 'pdom1');
model.result.dataset('dset3').set('solution', 'sol3');
model.result.dataset('cpt1').set('probetag', 'pdom1');
model.result.dataset('cpt1').set('data', 'dset3');
model.result.dataset('dset4').set('solution', 'none');
model.result.dataset('dset5').set('solution', 'none');
model.result.dataset('dset6').set('solution', 'none');
model.result.dataset('dset7').set('solution', 'none');
model.result.dataset('dset8').set('solution', 'none');
model.result.dataset('dset9').set('solution', 'sol3');
model.result.dataset('dset10').set('solution', 'none');
model.result.numerical.create('max1', 'MaxSurface');
model.result.numerical.create('min1', 'MinSurface');
model.result.numerical('max1').set('data', 'dset9');
model.result.numerical('max1').selection.set([1]);
model.result.numerical('max1').set('probetag', 'none');
model.result.numerical('min1').set('data', 'dset9');
model.result.numerical('min1').selection.set([1]);
model.result.numerical('min1').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result.create('pg3', 'PlotGroup2D');
model.result.create('pg4', 'PlotGroup2D');
model.result.create('pg5', 'PlotGroup2D');
model.result.create('pg6', 'PlotGroup2D');
model.result.create('pg7', 'PlotGroup2D');
model.result.create('pg8', 'PlotGroup2D');
model.result.create('pg10', 'PlotGroup1D');
model.result.create('pg11', 'PlotGroup2D');
model.result.create('pg12', 'PlotGroup2D');
model.result.create('pg13', 'PlotGroup2D');
model.result.create('pg14', 'PlotGroup2D');
model.result.create('pg15', 'PlotGroup2D');
model.result.create('pg16', 'PlotGroup2D');
model.result.create('pg17', 'PlotGroup2D');
model.result.create('pg18', 'PlotGroup2D');
model.result.create('pg19', 'PlotGroup1D');
model.result.create('pg20', 'PlotGroup1D');
model.result.create('pg21', 'PlotGroup2D');
model.result.create('pg22', 'PlotGroup2D');
model.result.create('pg23', 'PlotGroup2D');
model.result.create('pg24', 'PlotGroup2D');
model.result('pg1').set('data', 'dset9');
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
model.result('pg10').set('probetag', 'window1_default');
model.result('pg10').create('tblp1', 'Table');
model.result('pg10').feature('tblp1').set('probetag', 'pdom1/ppb1');
model.result('pg11').set('data', 'dset7');
model.result('pg11').create('surf1', 'Surface');
model.result('pg11').create('str1', 'Streamline');
model.result('pg12').set('data', 'dset7');
model.result('pg12').create('surf1', 'Surface');
model.result('pg12').create('str1', 'Streamline');
model.result('pg12').feature('surf1').set('expr', 'I');
model.result('pg13').set('data', 'dset7');
model.result('pg13').create('surf1', 'Surface');
model.result('pg13').create('str1', 'Streamline');
model.result('pg13').feature('surf1').set('expr', 'R');
model.result('pg14').set('data', 'dset7');
model.result('pg14').create('surf1', 'Surface');
model.result('pg14').create('str1', 'Streamline');
model.result('pg14').feature('surf1').set('expr', 'C');
model.result('pg15').set('data', 'dset9');
model.result('pg15').create('surf1', 'Surface');
model.result('pg15').create('str1', 'Streamline');
model.result('pg16').set('data', 'dset9');
model.result('pg16').create('surf1', 'Surface');
model.result('pg16').create('str1', 'Streamline');
model.result('pg16').feature('surf1').set('expr', 'I');
model.result('pg17').set('data', 'dset9');
model.result('pg17').create('surf1', 'Surface');
model.result('pg17').create('str1', 'Streamline');
model.result('pg17').feature('surf1').set('expr', 'R');
model.result('pg18').set('data', 'dset9');
model.result('pg18').create('surf1', 'Surface');
model.result('pg18').create('str1', 'Streamline');
model.result('pg18').feature('surf1').set('expr', 'C');
model.result('pg19').create('tblp1', 'Table');
model.result('pg20').create('tblp1', 'Table');
model.result('pg21').set('data', 'dset10');
model.result('pg21').create('surf1', 'Surface');
model.result('pg21').create('str1', 'Streamline');
model.result('pg22').set('data', 'dset10');
model.result('pg22').create('surf1', 'Surface');
model.result('pg22').create('str1', 'Streamline');
model.result('pg22').feature('surf1').set('expr', 'I');
model.result('pg23').set('data', 'dset10');
model.result('pg23').create('surf1', 'Surface');
model.result('pg23').create('str1', 'Streamline');
model.result('pg23').feature('surf1').set('expr', 'R');
model.result('pg24').set('data', 'dset10');
model.result('pg24').create('surf1', 'Surface');
model.result('pg24').create('str1', 'Streamline');
model.result('pg24').feature('surf1').set('expr', 'C');

model.component('comp1').probe('pdom1').genResult([]);

model.study('std2').label('Steady QS solver 1');
model.study('std2').feature('stat').set('useinitsol', true);
model.study('std2').feature('stat').set('initmethod', 'sol');
model.study('std2').feature('stat').set('initstudy', 'std2');
model.study('std2').feature('stat').set('solnum', 'last');
model.study('std2').feature('stat').set('usesol', true);
model.study('std2').feature('stat').set('notsolmethod', 'sol');
model.study('std3').label('Unsteady QS solver');
model.study('std3').feature('time').set('useinitsol', true);
model.study('std3').feature('time').set('initmethod', 'sol');
model.study('std3').feature('time').set('initstudy', 'std2');
model.study('std3').feature('time').set('solnum', 'last');
model.study('std3').feature('time').set('usesol', true);
model.study('std3').feature('time').set('notsolmethod', 'sol');

model.sol('sol2').attach('std2');
model.sol('sol2').feature('st1').label('Compile Equations: Stationary');
model.sol('sol2').feature('v1').label('Dependent Variables 1.1');
model.sol('sol2').feature('v1').set('control', 'user');
model.sol('sol2').feature('v1').set('initmethod', 'sol');
model.sol('sol2').feature('v1').set('initsol', 'sol3');
model.sol('sol2').feature('v1').set('solnum', 'last');
model.sol('sol2').feature('v1').set('notsolmethod', 'sol');
model.sol('sol2').feature('s1').label('Stationary Solver 1.1');
model.sol('sol2').feature('s1').feature('dDef').label('Direct 2');
model.sol('sol2').feature('s1').feature('dDef').set('thresh', 0.1);
model.sol('sol2').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol2').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol2').feature('s1').feature('fc1').set('initstep', 0.01);
model.sol('sol2').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol2').feature('s1').feature('fc1').set('maxiter', 50);
model.sol('sol2').feature('s1').feature('d1').label('Direct 1.1');
model.sol('sol2').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol2').feature('s1').feature('d1').set('pivotperturb', 1.0E-13);

model.sol('sol3').attach('std3');
model.sol('sol3').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol3').feature('v1').label('Dependent Variables 1.1');
model.sol('sol3').feature('v1').set('initmethod', 'sol');
model.study('std3').feature('time').set('initstudy', 'zero');
model.sol('sol3').feature('v1').set('notsolmethod', 'sol');

model.sol('sol3').feature('t1').label('Time-Dependent Solver 1.1');

model.sol('sol3').feature('t1').set('rtol', 0.005);
model.sol('sol3').feature('t1').set('maxorder', 2);
model.sol('sol3').feature('t1').set('stabcntrl', true);
model.sol('sol3').feature('t1').feature('dDef').label('Direct 2');
model.sol('sol3').feature('t1').feature('aDef').label('Advanced 1');
model.sol('sol3').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol3').feature('t1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol3').feature('t1').feature('fc1').set('linsolver', 'd1');
model.sol('sol3').feature('t1').feature('fc1').set('maxiter', 8);
model.sol('sol3').feature('t1').feature('fc1').set('damp', 0.9);
model.sol('sol3').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol3').feature('t1').feature('fc1').set('stabacc', 'aacc');
model.sol('sol3').feature('t1').feature('fc1').set('aaccdim', 5);
model.sol('sol3').feature('t1').feature('fc1').set('aaccmix', 0.9);
model.sol('sol3').feature('t1').feature('fc1').set('aaccdelay', 1);
model.sol('sol3').feature('t1').feature('d1').label('Direct, concentrations (tds) (merged) (merged) (merged)');
model.sol('sol3').feature('t1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol3').feature('t1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol3').feature('t1').feature('i1').label('AMG, concentrations (tds4)');
model.sol('sol3').feature('t1').feature('i1').set('maxlinit', 50);
model.sol('sol3').feature('t1').feature('i1').feature('ilDef').label('Incomplete LU 1');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').label('Multigrid 1.1');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').set('maxcoarsedof', 50000);
model.sol('sol3').feature('t1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol3').feature('t1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').label('SOR Line 1.1');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('iter', 1);
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.7);
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('relax', 0.5);
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').label('SOR Line 1.1');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('iter', 1);
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.7);
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('relax', 0.5);
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol3').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);

model.variable('var1').set('kp', 'kpb + kA *(x/L)');
model.sol('sol3').feature('v1').set('clist', {'range(0,1e6,1e6)' '150000.0[s]'});
model.sol('sol3').feature('t1').set('tlist', 'range(0,1e6,1e6)');
model.sol('sol3').feature('t1').set('maxstepconstraintbdf', 'auto');
disp('Performing initial time-dependent solve...')
model.sol('sol3').runAll;
disp('Performing initial steady solve...')
model.sol('sol2').runAll;
model.result.create('pg10', 'PlotGroup1D');
for T=Tsweep
    for K=Ksweep

model.param.set('T', num2str(T));
model.param.set('K', num2str(K));
model.sol('sol3').feature('v1').set('clist', {'range(0,T/100,T*10)' '150000.0[s]'});
model.sol('sol3').feature('t1').set('tlist', 'range(0,T/100,T*10)');
model.variable('var1').set('kp', 'kpb + kA * (x/L)*(1+K*sin(2*pi*t/T))');
model.sol('sol3').feature('v1').set('initsol', 'sol2');
model.sol('sol3').feature('v1').set('solnum', 'last');
model.sol('sol3').feature('t1').set('maxstepconstraintbdf', 'const');
model.sol('sol3').feature('t1').set('maxstepbdf', 'T/10');
disp('Performing full time-dependent solve...')
model.sol('sol3').runAll;

model.result.dataset('dset3').label('Probe Solution 3');
model.result.dataset('dset3').set('frametype', 'spatial');
model.result.dataset('cpt1').set('data', 'dset3');
model.result.dataset.remove('dset11');
model.result.numerical('max1').set('table', 'tbl1');
model.result.numerical('max1').set('expr', {'Q'});
model.result.numerical('max1').set('unit', {''});
model.result.numerical('max1').set('descr', {''});
model.result.numerical('min1').set('table', 'tbl6');
model.result.numerical('min1').set('expr', {'Q'});
model.result.numerical('min1').set('unit', {''});
model.result.numerical('min1').set('descr', {''});
model.result.numerical('max1').setResult;
model.result.numerical('min1').setResult;

Qmax = model.result.table('tbl1').getReal;
Qmin = model.result.table('tbl6').getReal;
save(['output_crit_T_',num2str(T),'_K_',num2str(K),'_kA_',num2str(kA),'.mat'],'Qmax','Qmin')
flag = 0;

if max(Qmax(401:end,2))>1
flag = 1;
    break
end
    
    end
if flag==1
break
end
end

Tc = T;
