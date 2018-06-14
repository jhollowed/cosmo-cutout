# cosmo-cutout

Modeling synthetic observations of the sky is an important task in modern cosmology, with the general hierarchy of necessary data products being as follows:

* Develop and run simulations of structure formation within a cosmological volume (from this we obtain the mass distribution of the contents of the universe). 
* Place an "obeserver" at some location within the simulation volume, and construct a *lightcone*, which is a four dimensional hypersurface on which lie all spacetime events output by the simulation that are seen by the observer at *z*=0 (we then have the mass distribution that can actually be seen in the sky today). 
* Finally, extract a particular field fo view from this lightcone that a telescope would actually observe during an exposure.
 
The code distributed in this repository performs the final step above, slicing out a rectangular cone, from a cosmological lightcone catalog, defined by some angular bounds with respect to an observer.

## Getting Started

### Prerequisites

This project utilizes [GenericIO](https://trac.alcf.anl.gov/projects/genericio) to read in simulation lightcone outputs, as well as [OpenMP](https://www.openmp.org/) and [MPI](https://www.mpi-forum.org/). The code is intended to be run on the sytem environments at the [ALCF](https://www.alcf.anl.gov/computing-resources), and will otherwise need its Makefile to be edited; specifically, look toward the `TRUNK`, `LIBS` and `INCLUDES` to ensure that GenericIO is properly linked. The next section of this readme will assume that [HACC](https://arxiv.org/abs/1410.2805), or some other source of GenericIO, is already built. If using an ALCF machine, the `Makefile`s present in the repo already point to the correct locations.

For help on running at ALCF, there are pages dedicated to compiling and linking on [BG/Q systems](https://www.alcf.anl.gov/user-guides/bgq-systems-compiling-linking) and [Cooley](https://www.alcf.anl.gov/user-guides/cooley-compling-linking).

### Installing

If running at ALCF, first clone the repo

```
git clone https://github.com/jhollowed/cosmo-cutout.git
```

and make the following environment changes...
First edit your `~/.soft`

on Mira/Cetus:
```
+mpiwrapper-xl.legacy.ndebug
@default
```

on Cooley:
```
+mvapich2
@default
```

and run 

```
resoft
```

Or, if running on Pheonix (formerly Jupiter), run

```
module load mvapich2-2.2b-gcc-5.3.0-o4of6w7
```

Next, compile the source code

on Mira/Cetus:
```
make -f Makefile.mira
```

on Cooley:
```
make -f Makefile.cooley
```

on Pheonix:
```
make -f Makefile.pheonix
```

## Running the cutout

*Notation:*

&#x03B8; - coaltitutde coordinate

&#x03D5; - azimuthal coordinate

*x*, *y*, *z* - comoving cartesian coordinates

`lc_cutout` - the executable resultant from following the installation instructions

Two use cases are supported: 

### Use case 1: Constant angular bounds

Define the theta and phi bounds explicitly:
  
```
lc_cutout <input lightcone directory> <output directory> <depth> --theta <theta_center> <d_theta> --phi <phi_center> <d_phi;>
```

where the &#x03D5;<sub>center</sub> argument is the azimuthal coordinate of the center of the field of view that one wishes to cut out of the lightcone, and d&#x03D5; is the angualar distance from this center to the edge of the cutout, and likewise for the similar &#x03B8;  args. That is, the result will be a sky area that spans 
 
(&#x03B8;<sub>center</sub> - d&#x03B8;) < &#x03B8; < (&#x03B8;<sub>center</sub> + d&#x03B8;)

(&#x03D5;<sub>center</sub> - d&#x03D5;) < &#x03D5; < (&#x03D5;<sub>center</sub> + &#x03D5;)

The expected angular units are DEGREES. The `--theta` and `--phi` flags can be replaced with `-t` and `-p`.

The `depth` parameter is the maximum desired redshift of the cutout (limited, of course, by the maximu redshift of the input lightcone catalog).


### Use case 2: Nonlinear angular bounds

Allow the &#x03B8; and &#x03D5; bounds to be computed interanally to obtain a cutout of a certain width (*box length*), in Mpc/h, centered on a certain cartesian positon, (*x*&#x2080;, *y*&#x2080;, *z*&#x2080;) Mpc/h (intended to be used for making cutouts centerd on specific simulation objects, like halos):

```
lc_cutout <input lightcone directory> <output directory> <depth> --halo <x_0> <y_0> <z_0> --boxLength <box length>
```

where the `--halo` and `--boxLength` flags can be replaced with `-h` and `-b`. Those options give the position of the object to center the cutout on, and the one-dimensional size of the field of view in Mpc/h, tangent to the line of sight, at the distance of the object of interest, respectively.

We want to express the positions of  all of our lightcone objects in spherical coordinates, to perform the cutout, and we want that coordinate system to be rotated such that the object of intererst lies on the equator at

(*r*, 90&#x00B0;, 0&#x00B0;)

where 

*r* = (*x*<sup>2</sup> + *y*<sup>2</sup> + *z*<sup>2</sup>)<sup>1/2</sup>;

Let's call the position vector of the halo before this rotation 

**a** = [*x*&#x2080;, *y*&#x2080;, *z*&#x2080;], 

and after,

 **b** = [*x*<sub>rot</sub>, *y*<sub>rot</sub>, *z*<sub>rot</sub>] = [*r*, 0, 0]

We perform this rotation for each lightcone object via the *Rodrigues rotation formula*, which answers the following question: given a position vector **v**, a normalized axis of rotation **k**, and an angle of rotation &#x03B2;, what is an analytical form for a new vector **v**<sub>rot</sub> which is **v** rotated by an anlge &#x03B2; about **k**?

First, we find **k** by taking the cross product of two vectors defining the 
plane of rotation. The obvious choice of these two vectors are **a** and **b**, as 
defined above;

k = (**a** &#x2A2F; **b**) / &#x2016;**a** &#x2A2F; **b**&#x2016;

then, for any other position vector **v**, **v**<sub>rot</sub> is given by

**v**<sub>rot</sub> = **v**cos&#x03B2; + (**k** &#x2A2F; **v**)sin&#x03B2; + **k**(**k** &#x22C5; **v**)(1-cos&#x03B2;)

This coordinate rotation is required because the bounding functions which define the field of view of the observer, while constant in theta-phi space, are nonlinear in cartesian space. The distortion is maximized near the poles of the spherical coordinate system, and minmized at the equator. Areas defined by constant theta-phi bounds then appear trapezoidal to the observer when far from the equator. It is important that our cutout areas are maintained as square for at least two reasons:

* At the moment, FFT restrictions of flat-sky lensing code require that the cutout is square
* The cutouts returned will not actually have all side lengths of `boxLength` if we don't do this rotation, which the user explicitly requested

Note that the first use case describe does not perform the coordinate rotation which is described in the second. So, cutouts returned will not necessarily be square, or symmetrical.

## Authors

wip

## Acknowledgments

wip

