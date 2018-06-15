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

# Running the cutout

### *Notation:*

*&#x03B8;* - coaltitutde coordinate

*&#x03D5;* - azimuthal coordinate

*x*, *y*, *z* - comoving cartesian coordinates

`lc_cutout` - the executable resultant from following the installation instructions

Two use cases are supported: 

## Use case 1: Constant angular bounds

Define the theta and phi bounds explicitly:
  
```
lc_cutout <input lightcone directory> <output directory> <min redshift> <max redshift> --theta <theta_center> <d_theta> --phi <phi_center> <d_phi;>
```

#### Arguments:

  * `input lightcone directory` - the location of the top-level directory of a simulation lightcone. The directory structure is expected to match that as described in section 4.5 (Figure 7) of [Creating Lightcones in HACC](http://www.joehollowed.com/lightcone_notes.html), which contains step-wise subdirectories  
  * `output directory` - where to save the result of the cutout. A new subdirectory will be created at this location, of the form `lcCutoutXXX` for each step `XXX` included in the calculation  
  * `min redshift` - allows the user to begin the cutout construction at some distance away from the observer position (the origin) in redshift-space. Setting this parameter to a nonzero value is intended to be used in the case that the user is breaking up the cutout procedure across separate jobs  
  * `max redshift` - controls the depth of the cutout in redshift-space (limited, of course, by the maximum redshift of the input lightcone catalog)  
  * `theta_center` - the &#x03B3; coordinate of the center of the desired cutout field of view  
  * `d_theta` - the distance from &#x03B3;<sub>center</sub> to the edge of the field of view  
  * `phi_center` the &#x03D5; coordinate of the center of the desired cutout field of view  
  * `d_phi` - the distance from &#x03D5;<sub>center</sub> to the edge of the field of view  
 
That is, the result of the cutout will be a sky area that spans 
 
(&#x03B8;<sub>center</sub> - d&#x03B8;) < &#x03B8; < (&#x03B8;<sub>center</sub> + d&#x03B8;)

(&#x03D5;<sub>center</sub> - d&#x03D5;) < &#x03D5; < (&#x03D5;<sub>center</sub> + d&#x03D5;)

and redshift range...

The expected angular units are DEGREES. The `--theta` and `--phi` flags can be replaced with `-t` and `-p`.

## Use case 2: Nonlinear angular bounds

Allow the &#x03B8; and &#x03D5; bounds to be computed interanally to obtain a cutout of a certain width (*box length*), in Mpc/h, centered on a certain cartesian positon, (*x*&#x2080;, *y*&#x2080;, *z*&#x2080;) Mpc/h (intended to be used for making cutouts around specific simulation objects, like halos):

```
lc_cutout <input lightcone directory> <output directory> <min redshift> <max redshift> --halo <x_0> <y_0> <z_0> --boxLength <box length>
```

#### Arguments:

  * `input lightcone directory`, `output directory`, `min redshift`, `max redshift` - See description above  
  * `x_0`, `y_0`, `z_0` - The comoving cartesian position, in Mpc/h, of the object on which to center the cutout  
  * `box length` - the width of the fov around the object of iterest, in Mpc/h at the distance of the object (let this value be denoted as *B*, then *d&#x03B8*, as defined above, is tan<sup>-1</sup>(*B*/2*r*), where *r* is *r* = (*x*<sub>0</sub><sup>2</sup> + *y*<sub>0</sub><sup>2</sup> + *z*<sub>0</sub><sup>2</sup>)<sup>1/2</sup>)  

The `--halo` and `--boxLength` flags can be replaced with `-h` and `-b`.

### Multiple objects of interest

If one has many objects of interest around which they would like lightcone cutouts, then it would be inefficient to call the above command wiht the `-h` option each time, since each one of those runs would be using resources to re-read the same input lightcone data (which has the potential to be very large). To address this, the program can be run in the following manner:

```
lc_cutout <input lightcone directory> <output directory> <min redshift> <max redshift> --haloFile <input object file> --boxLength <box length>
```

#### Arguments:

  * `input lightcone directory`, `output directory`, `min redshift`, `max redshift`, `box length` - See description above.  
  * `input object file` - A plain text file containing one line per object of interest, which includes an identifying object tag, and three cartesian comoving positions, as such:  
  
```
123456789 50 55 20
987654321 10 0 30
...
```
where the object tags are expected to be of type `int64_t` (satisfied by HACC particle id's and halo tags). In this example, the first object has a tag of `123456789`, and position *x*=`50`, *y*=`55`, *z*=`20`.

The `--haloFile` option can also be specified with `-f` (and, as above, `--boxLength` with `--b`). The cutouts for each of these objects will now be performed serially, with the lightcone read-in happening only once.


<details><summary>Click here to expand details on how exactly the cutout computation is done for Use Case 2</summary>
<p>

We want to express the positions of  all of our lightcone objects in spherical coordinates, to perform the cutout, and we want that coordinate system to be rotated such that the object of intererst (which we will from now on assume is a halo) lies on the equator at

(*r*, 90&#x00B0;, 0&#x00B0;)

where 

*r* = (*x*<sub>0</sub><sup>2</sup> + *y*<sub>0</sub><sup>2</sup> + *z*<sub>0</sub><sup>2</sup>)<sup>1/2</sup>

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
</p>
</details>

## Caveats

* Note that the Use Case 1 does not perform the coordinate rotation which is described in Use Case 2 (under the "click to expand" details). So, cutouts returned will not necessarily be square, or symmetrical, is far from the coordinate equator.

* The parallelism in this application occurs *spatially*, not temporally. That is, the lightcone volume is decomposed across *MPI* ranks, which prallelizes the read in, computation, and write-out. But there is no parallelism in *redshift*-space, meaning that each lightcone "step" (portion of the lightcone volume originating from a particular simulation snapshot) are treated in serial. Further, if option `-f` is used as described under Use Case 2, then those multiple requested cutouts are also treated serially. 

* The requested `min redshift` and `max redshift` are converted to a simulation step number assuming a simulation run that included 500 total time steps, and began at a redshift of 200. At the moment, there is no way for the user to easily change this, other than modifying the calls to `getLCSteps()` in `src/main.cpp` and rebuilding.


#  Authors

wip

## Acknowledgments

wip

