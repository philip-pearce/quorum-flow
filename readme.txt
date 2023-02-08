For paper "Emergent robustness of bacterial quorum sensing in fluid flow"

Comsol simulations are performed via the files solve_steady_comsol.m and solve_steady_comsol_3D.m.
This requires Comsol Livelink with Matlab (Comsol version 5.5).

To access a Comsol simulation in the graphical user interface, run one of these simulations in Matlab and then
type

model=ModelUtil.model('Model')

on the Matlab command line. Then the model can be saved by typing

model.save('insert_path_here')

To recreate the 2D results, run write_parameters.m to generate parameter sets,
then sweep_density_Pe_params.m to generate results. First create a folder
called 'res' in which to save the results.

To recreate the 3D results, run sweep_density_3D.m. In these simulations,
Comsol files are saved to disk for further analysis (they are quite large).
First create a folder called 'res3D' in which to save the results.

To recreate the results for the time-dependent figures, open clean_td.mph and
run the time-dependent simulations with the parametric sweep defining the values
of the time period 1/epsilon. (Note that first a steady solve should be performed
for the mean Peclet number; when performing this steady solve, the variable Pe must
be edited to be equal to Pe_mean; then add the time-dependent part back in for the
unsteady solve).

For paper "Universal dynamics of biological pattern formation in spatio-temporal morphogen variations"

Comsol simulations performed in Comsol Livelink with Matlab (Comsol version 5.5) as above.
To run the simulations, specify the oscillation time period (T_sweep), the required value(s) of K [oscillation size;
see Fig. S1 in paper] and the kinetic parameter kA [average value of parameter that is oscillating in time].
The code will output a Matlab file with the maximum and minimum concentrations of the signalling molecule against time.
The code will stop for the value of T at which the maximum concentration reaches a critical value
(see Fig. S1 in paper).
