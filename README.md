saxs_cylinder_simulation
========================

The code for programs used in obtaining the simulation results in the paper H. Suhonen et al. Phys. Med. Biol 50 (2005), 5401-5416. 


Essentially the program calculates small angle scattering from a set of
cylinders that have axially periodic electron densities. This package
consists of two programs: SAXS_calc calculates the scattering patterns
and cylgen generates the sets of cylinders that SAXS_calc uses.



SAXS_calc
=========

Usage
-----
SAXS_calc -k kValueFile [-c axialDensityFile] scatteringObjectFile



Arguments
---------
SAXS_calc is the program that calculates the scattering. The input
options are the following


	-k kValueFile

	This gives the name of a file where the k-values at which to
	calculate the scattering are listed. The k-value file should
	consists of triplets of numbers (each triplet on its own
	line). The values give the x-,y- and z-components of the
	scattering vector in the units of (1/Å). x-direction is the
	direction from which the beam is coming from and y- and
	z-vectors are in the detector plane. The z-direction is the
	meridional direction for a cylinder whose polar angle is
	0. The definition of k that is used is k = 4*pi*sin(theta) /
	lambda. 


	-c axialDensityFile

	This can be used to define the axial electron density of the
	cylinder more accurately than with just a step function. The
	file consists of lines with two numbers. The first number
	gives the relative distance along the period (between 0 and 1)
	and the second number gives the electron density value at that
	point. 



Input
-----

scatteringObjectFile should be in the format produced by the program
cyl_gen. This file contains a description of the cylinders that are
used for calculation.


Output
------
SAXS_calc produces output in the format:

<pre>
unit_count = 1000.000000
cylinders_per_unit_min = 1.000000
cylinders_per_unit_max = 7.000000
r_avg = 300.000000
r_stdev = 30.000000
l_avg = 650.000000
l_stdev = 3.000000
alpha_stdev = 10.000000
zdisp = 0.000000
rho1 = 0.500000
rho2 = 0.700000
period_proportion = 0.500000
period_proportion_stdev = 0.000000
period_min = 7.000000
period_max = 17.000000
tilt_lattice = 1.000000
use_random_orientation = 1.000000
k file name = kfile_eq
electron density file name =
*********
0.000000 0.000000 0.000000 36472725433394957123584.0000000000000000
0.000000 0.000100 0.000000 35598907331799645421568.0000000000000000
0.000000 0.000200 0.000000 33157384857680483450880.0000000000000000
0.000000 0.000300 0.000000 29625397241660145926144.0000000000000000
0.000000 0.000400 0.000000 25618389524108539854848.0000000000000000
0.000000 0.000500 0.000000 21704625012988807479296.0000000000000000
0.000000 0.000600 0.000000 18268495742396761374720.0000000000000000
</pre>

First is a header part that lists the parameters that were used.
After the header comes the actual scattering data. A row consists of
the three components of the k-vector values (4*pi*sin(theta)/lambda)
followed by the intensity value:

k_x k_y k_z intensity

The intensity values are not normalized (and the absolute magnitude of
these values does not relate to any property of real collagen systems,
but reflects only the number and total volume of the cylinders used).




cylgen
=======


Usage
-----
cylgen description_file



Arguments
---------

none




Input
-----

A file that describes the properties for the cylinders to be
created. The following parameters can be used in this file:

unit_count
the number of independently scattering cylinder units


cylinders_per_unit_min
cylinders_per_unit_max
The limits for number of cylinders in a scattering unit. The number of
cylinders for each scattering unit will be drawn uniformly from this
interval. Interference effects within each scattering unit will be
included in the calculation, but interference effects between
different units are ignored.

r_avg
r_stdev
The mean and standard deviation for the gaussian radius distribution
of the cylinders (in units of Å). 

l_avg
l_stdev
The mean and standard deviation for the gaussian length distribution
of the cylinders (i.e. the length of single period) (in units of Å). 


alpha_stdev
The width of the polar angle distribution (standard deviation of a
gaussian distribution) in units of degrees.


zdisp
Number of Å the cylinders within a bunch could be displaced axially
relative to each other (recommended value = 0).

rho1
rho2
The relative electron densities of the two parts in the step-function
axial electron density.


period_proportion
period_proportion_stdev
The mean relative length of the first part (relative to total) in the step-function axial
electron density and its standard deviation (gaussian distribution).

period_min 
period_max 
Minimun and maximum number of periods in a cylinder. The number of
periods for each cylinder separately will be selected uniformly from
this range.


tilt_lattice
Tilt also the lattice when tilting the cylinders (0 or 1). Recommended
value is 1.


use_random_orientation
Value 0 or 1. If 1 is used, alpha_stdev has no effect and the cylinder
orientation is selected randomly.


packdist_radius_fraction
packdist_addition
These control how closely the cylinders are packed to each
other. Values of 1 and 0 (respectively) give the closest possible
packing. Increasing packdist_radius_fraction increases the packing
distance by a fraction of the cylinder radii. Increasing
packdist_addition adds a uniformly distributed random distance to the
packing. To see how these parameters (and the packing) work, it may be
best to visualize the cylinders and try with different parameter
values. 



Output
------

	     
The output is generated on the standard output, which should be
redirected to the appropriate file. 


The output has a header part that includes the parameters that were
use in generating the data. After that follow cylinder descriptions
and each scattering unit (i.e. independent scatterers) is separated by
a line that contains a *.

<pre>
###
unit_count = 1000.000000
cylinders_per_unit_min = 1.000000
cylinders_per_unit_max = 7.000000
r_avg = 300.000000
r_stdev = 30.000000
l_avg = 650.000000
l_stdev = 3.000000
alpha_stdev = 10.000000
zdisp = 0.000000
rho1 = 0.500000
rho2 = 0.700000
period_proportion = 0.500000
period_proportion_stdev = 0.000000
period_min = 7.000000
period_max = 17.000000
tilt_lattice = 1.000000
use_random_orientation = 1.000000
###
cylper 0.099435 0.010617 0.000000 250.929376 69.983011 276.094774 651.564242 14 0.500000 0.700000 0.500000
cylper -522.928788 -55.836727 0.000000 275.072001 69.983011 276.094774 647.163534 12 0.500000 0.700000 0.500000
cylper -253.624547 157.301014 503.258329 334.207423 69.983011 276.094774 652.049472 17 0.500000 0.700000 0.500000
*
cylper 0.009597 0.099538 0.000000 301.215310 48.552715 354.492964 653.681838 14 0.500000 0.700000 0.500000
*
cylper -0.096578 0.025935 0.000000 331.523882 39.929811 74.968435 650.226145 9 0.500000 0.700000 0.500000
cylper 596.851067 -160.278230 0.000000 286.573130 39.929811 74.968435 643.701919 10 0.500000 0.700000 0.500000
cylper 445.844223 297.682201 -337.422733 301.971621 39.929811 74.968435 647.048766 17 0.500000 0.700000 0.500000
cylper -76.367907 450.174688 -347.331604 242.140182 39.929811 74.968435 651.313901 17 0.500000 0.700000 0.500000
cylper -550.447423 296.900070 -120.514898 305.299598 39.929811 74.968435 649.365565 10 0.500000 0.700000 0.500000
cylper 245.477975 -433.477379 297.123416 248.574652 39.929811 74.968435 651.080265 11 0.500000 0.700000 0.500000
*
</pre>

The parameters for each cylinder are in the following order:

x y z r alpha beta l n rho1 rho2 period_proportion

x,y,z = the coordinates of the cylinder center point (in Å)

r = the radius of the cylinder

alpha = the polar angle with respect to the z-axis (degrees)

beta = azimuthal angle (i.e. rotation around z) measured against
x-axis (degrees)

l = the length of a period

n = the number of periods

rho1 = electron density 1

rho2 = electron density 2

period_prortion = length of the rho1 region as a fraction of the total
period length
