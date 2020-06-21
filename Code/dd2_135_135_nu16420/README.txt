The format of the data is guided by previous outputs of NS-NS merger
models of Bauswein et al.
(see https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/Just2015/index.html) with some changes marked with "*".

1- The trajectory.datXXXXXX files provide data for a given particle trajectory. 
The suffix is the number of the ejected SPH particle.

The first line is a header with the following columns:

 - the number (ID) of the SPH particle
 - the mass of the SPH particle in solar masses
 - the initial electron fraction*
 - the time until when the trajectory evolution is reliable (end of simulation) in ms
 - the initial rest-mass density in g/cm^3
 - the maximum temperature of the trajectory during its evolution in MeV
 - dummy*

The following lines provide the time evolution of trajectory 
with the following entries in columns:

    - Time in ms
    - Rest-mass density in g/cm^3
    - Pressure in dyne/cm^2
    - Electron fraction*
    - Velocity of the SPH particle in km/s
    - Temperature in MeV (not postprocessed!)*
    - Electron chemical potential (including electron rest mass)*
    - Proton chemical potential (without rest mass)*
    - Neutron chemical potential (without rest mass)*
    - Radial coordinate of the particle (in Km)*
    - Angular coordinate of the particle (in deg)*
The Fortran code used to write the files is

write(908,fmt='(I6,6e12.4)')particle_ID,particle_mass,ye0,end_time,rho0,T_max,theta_final

for the header and 

write(908,fmt='(11e12.4)')time(k), rho(k), P(k), ye(k), v(k),
temp(k),mue(k),mup(k),mun(k), r(k),theta(k)

for the body.

2 - Now the .tar files also contain an additional file named
"neutrino_properties_bins.dat" which harbours the neutrino luminosities
and mean energies in angular bins, measured at 50 and 100 Km, as a
function of time. The first line in the file contains two integers with
the number of angular bins (18) and the angular width of each bin (10
deg). The rest of the file is written as:

open(unit=3453,file='neutrino_properties_bins.dat',status='unknown',position='append',form='formatted',recl=1024)
write(3453,*)time,lumibin50(1:nanglebins,1) , lumibin50(1:nanglebins,2) ,
lumibin100(1:nanglebins,1), lumibin100(1:nanglebins,2), 
meanebin50(1:nanglebins,1), meanebin50(1:nanglebins,2),
meanebin100(1:nanglebins,1), meanebin100(1:nanglebins,2)
close(3453)

Where time (real) is in ms (measured from the beginning of the simulation),
"lumibin" (double precision) contains the neutrino energy luminosities
(erg/s) and "meanebin" (double precision) contains the mean neutrino
energies (MeV). There are two sets of arrays containing angular bins,
distinguished by the values 50 or 100 in the variable names. Each of these
sets contains the neutrino properties measured at 50 or 100 km from the
center of the remnant. Each array has then two indexes: the bin number,
and an integer value, 1 for nue and 2 for nuebar.

Note: the time (in ms) is not normalized in any of the files. The first 
entry (timestep) in each of he files corresponds to the time of minimum 
lapse, and thus the 0 in the plots in Ardevol-Pulpillo et al (2019).

Please, do not hesitate to contact Ricard Ardevol-Pulpillo 
(ricardar@mpa-garching.mpg.de) if there is any issue with the files
or the formatting.


Details concerning the physics used in the models and the data extraction:
--------------------------------------------------------------------------

First of all, the simulations were performed exactly as described in
Ardevol-Pulpillo et al (2019). The Cartesian grid covering the NSs in
which the neutrino interactions are calculated expands 100km in all three
directions, with a resolution of 0.5km. After the merger, as the remnant
torus expands, the grid size is progressively increased keeping the whole
remnant covered at all times, keeping always the same resolution (cell
size) throughout the simulation.

The extraction of the ejecta properties was performed exactly in the same
way as previous publications, with the exception that we chose NOT to post
process the temperature. The reasoning is simple: we now evolve the NS
matter including weak interactions which drive the evolution of the
temperature (through the internal energy) and are highly sensitive to it.
The post processing of the temperature would thus create an inconsistency
between the neutrino properties and the other thermodynamical properties.

The extraction of the neutrino properties contained in
"neutrino_properties_bins.dat" is done in a slightly different fashion to
the procedure detailed in Ardevol-Pulpillo et al (2019). This is namely
due to different requirements: in the paper we needed global values for
such quantities, while now we are interested in local values (angular and
radial dependency).

I would recommend the reading of section 2.3 from the paper to familiarize 
oneself with the ray tracing approach employed for the neutrino
reabsorption, which is also used to extract the neutrino properties. To
summarize in a few sentences, we assume neutrinos to be emitted as rays in
the direction of the gradient of the neutrino energy density. These rays
then travel through a Cartesian grid and deposit energy and lepton number
in the matter contained in the hit cells. We extract the neutrino
properties from these rays at the chosen radii (50 and 100 km) and then
store them in the corresponding angular bins.

The luminosity (L) emitted by a given cell is described in equation (34) as:

L=Q*V

with Q the neutrino energy source term and V the cell volume. This
luminosity is carried along the ray, redshifted and reduced due to absorption as
described in detail in section 2.3. When the ray crosses the spherical
surface at the chosen radius (50 or 100 km), we store the local value of
the ray's luminosity in the corresponding angular bin. So after all rays
are accounted for, the luminosity of a given angular bin is simply defined
by:

L_nu(r, theta) = sum( L(ray,r,theta) )

The neutrino mean energies are extracted in a similar spirit as the post
processed mean energies in section 2.5 from Ardevol-Pulpillo et al (2019):
a combination of optically thick and optically thin treatments.
- For rays originating or traversing optically thick conditions (optical
depth >2/3), the mean neutrino energies correspond to the local mean
energy of the neutrino spectrum:

E_nu^thick=T*ffen3/ffen2

where T is the femperature and ffenX the Fermi integral of order X. If the
ray leaves the optically thick regime it carries the neutrinospheric value
of the mean neutrino energies along its path, redshifted accordingly to
the corresponding radius (50 or 100 km) as described in Ardevol-Pulpillo
et al (2019), eq (48).

-For rays originating in the optically thin regime, we calculate the mean
energy in the leakage fashion:

E_nu^thin= Q/R

with R as the neutrino number source term. The ray then carries this mean
energy, redshifted accordingly, from the source to the corresponding
radius (50 or 100km). Note that if one such rays crosses a region with
optically thick conditions, the rules for optically thick conditions above
apply, as neutrinos will then thermalize.

The mean neutrino energy of an angle bin at a given radius (50 or 100) is
then calculated as the sum of the mean neutrino energies of the rays
crossing the angular surface of the bin, weighted by the neutrino energy
luminosity of each such rays.

E_nu(r, theta) = sum( E_nu(ray,r,theta) * L(ray,r,theta)) / L(r, theta)

Just as a final note for clarification, all luminosities and mean neutrino 
energies have been redshifted for a locar observer sititng at the 
corresponding surface at radius 50/100 km. Please bear in mind that the 
plots in Ardevol-Pulpillo et al (2019) such quantities were redshifted 
for an observer at infinity.

This covers all the new ingredients included in the trajectories.

Please do not hesitate to ask Ricard Ardevol-Pulpillo 
(ricardar@mpa-garching.mpg.de) or Thomas Janka (thj@mpa-garching.mpg.de)
if you have any questions or requests.





D
B
A
A

