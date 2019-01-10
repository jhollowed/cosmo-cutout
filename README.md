# cosmo-cutout

Modeling synthetic observations of the sky is an important task in modern cosmology, with the general hierarchy of necessary data products being as follows:

* Develop and run simulations of structure formation within a cosmological volume (from this we obtain the mass distribution of the contents of the universe). 
* Place an "observer" at some location within the simulation volume, and construct a *lightcone*, which is a four dimensional hypersurface on which lie all spacetime events output by the simulation that are seen by the observer at *z*=0 (we then have the mass distribution that can actually be seen in the sky today). 
* Finally, extract a particular field of view from this lightcone that a telescope would actually observe during an exposure.
 
The code distributed in this repository performs the final step above, slicing out a rectangular cone, from a cosmological lightcone catalog, defined by some angular bounds with respect to an observer.

## Getting Started

### Prerequisites

This project utilizes [GenericIO](https://trac.alcf.anl.gov/projects/genericio) to read in simulation lightcone outputs, as well as [OpenMP](https://www.openmp.org/) and [MPI](https://www.mpi-forum.org/). The code is intended to be run on the system environments at the [ALCF](https://www.alcf.anl.gov/computing-resources), and will otherwise need its Makefile to be edited; specifically, look toward the `TRUNK`, `LIBS` and `INCLUDES` to ensure that GenericIO is properly linked. The next section of this readme will assume that [HACC](https://arxiv.org/abs/1410.2805), or some other source of GenericIO, is already built. If using an ALCF machine, the `Makefile`s present in the repo already point to the correct locations.

For help on running at ALCF, there are pages dedicated to compiling and linking on [BG/Q systems](https://www.alcf.anl.gov/user-guides/bgq-systems-compiling-linking) and [Cooley](https://www.alcf.anl.gov/user-guides/cooley-compling-linking).

### Installing

If running at ALCF, first clone the repo

```bash
git clone https://github.com/jhollowed/cosmo-cutout.git
```

and make the following environment changes...

First edit your `~/.soft` to include these lines

on Mira/Cetus:
```bash
+mpiwrapper-xl.legacy.ndebug
@default
```

on Cooley:
```bash
+mvapich2-2.1-intel
+intel-composer-xe-2015
@default
```

and run 

```bash
resoft
```

Next, compile the source code

on Mira/Cetus:
```bash
make -f Makefile.mira
```

on Cooley:
```bash
make -f Makefile.cooley
```

# Running the cutout

### *Notation:*

*&#x03B8;* - coaltitutde coordinate

*&#x03D5;* - azimuthal coordinate

*x*, *y*, *z* - comoving Cartesian coordinates

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
  * `theta_center` - the *&#x03B8;* coordinate of the center of the desired cutout field of view  
  * `d_theta` - the distance from *&#x03B8;*<sub>center</sub> to the edge of the field of view  
  * `phi_center` the *&#x03D5;* coordinate of the center of the desired cutout field of view  
  * `d_phi` - the distance from *&#x03D5;*<sub>center</sub> to the edge of the field of view  
 
That is, the result of the cutout will be a sky area that spans 
 
(*&#x03B8;*<sub>center</sub> - *d&#x03B8;*) < *&#x03B8;* < (*&#x03B8;*<sub>center</sub> + *d&#x03B8;*)

(*&#x03D5;*<sub>center</sub> - *d&#x03D5;*) < *&#x03D5;* < (*&#x03D5;*<sub>center</sub> + *d&#x03D5;*)

and redshift range...

The expected angular units are DEGREES. The `--theta` and `--phi` flags can be replaced with `-t` and `-p`.

## Use case 2: Nonlinear angular bounds

Allow the *&#x03B8;* and *&#x03D5;* bounds to be computed internally to obtain a cutout of a certain width (*box length*), in arcminutes, centered on a certain Cartesian position, (*x&#x2080;*, *y&#x2080;*, *z&#x2080;*) Mpc/h (intended to be used for making cutouts around specific simulation objects, like halos):

```
lc_cutout <input lightcone directory> <output directory> <min redshift> <max redshift> --halo <x_0> <y_0> <z_0> --boxLength <box length>
```

#### Arguments:

  * `input lightcone directory`, `output directory`, `min redshift`, `max redshift` - See description above  
  * `x_0`, `y_0`, `z_0` - The comoving Cartesian position, in Mpc/h, of the object on which to center the cutout  
  * `box length` - the angular width of the fov around the object of interest, in arcmin (let this value be denoted as *B*, then *d&#x03B8;*, as defined above, *B*/2)  

The `--halo` and `--boxLength` flags can be replaced with `-h` and `-b`. 

### Multiple objects of interest

If one has many objects of interest around which they would like lightcone cutouts, then it would be inefficient to call the above command with the `-h` option each time, since each one of those runs would be using resources to re-read the same input lightcone data (which has the potential to be very large). To address this, the program can be run in the following manner:

```
lc_cutout <input lightcone directory> <output directory> <min redshift> <max redshift> --haloFile <input object file> --boxLength <box length>
```

#### Arguments:

  * `input lightcone directory`, `output directory`, `min redshift`, `max redshift`, `box length` - See description above.  
  * `input object file` - A plain text file containing one line per object of interest, which includes an object identifier, a halo redshfit, the snapshot/lightcone shell of that halo, mass, and optionally a radius, and three Cartesian comoving positions, as such:  
  
```
123456789 0.5, 1e14, 0.9, 50.0 55.0 20.0
987654321 0.1, 1e15, 1.4, 10.0 0.0 30.0
192837465 1, 1e13.5, 0.7, 110.0 35.0 20.0
...
```
In this example, the first object has an id of `123456789`, a redsshift of `0.5`, SO mass of `1e14 M_sun/h`, r*&#x2082;**&#x2080;**&#x2080;* radius of `0.9 Mpc/h`, position *x*=`50`, *y*=`55`, *z*=`20 Mpc/h`. All quantities expect for the `id` must be in such a form that they can be parsed as `floats`. The radius is optionally included, and can be removed as long as `massDef` is set to `fof` (see below). 

The identifiers can be anything, and are parsed as strings (in this way, they can be used for storing other meta data if desired). Under this usage, a new subdirectory will be created per object as listed in the `input object file` under `output directory`, of the form `halo_123456789`, for example. It is then under that directory that simulation step-wise directories will be created (as in the description of the `output directory` argument). It is also within the `output directory` that a `properties.csv` file will be written, which will contain the halo redshift, mass, optionally radius, and information about the scale of the final cutout. Any invalid/missing quantities in that file will be recoded as `-1` (this occurs when running Use Case 2 with `-h` rather than `-f`).

The `--haloFile` option can also be specified with `-f` (and, as above, `--boxLength` with `-b`). The cutouts for each of these objects will now be performed serially, with the lightcone read-in happening only once.

<details><summary>
<b><i>Click here to expand details on how exactly the cutout computation is done for Use Case 2</i></b>
</summary>
<p>

In making general lightcone cutouts around specific objects, a coordinate rotation is required. This is because the bounding functions which define the field of view of the observer, if they are described in constant angular terms, are in general nonlinear in Cartesian space (meaning the field of view will not appear square from the observer position). This distortion is maximized near the poles of the spherical coordinate system, and minimized at the equator, as long as the small-angle approximation holds for the size of the cutout. Areas defined by constant  *&#x03B8;* - *&#x03D5;* bounds, then, appear trapezoidal to the observer when far from the equator. It is important that our cutout areas are maintained as square for at least two reasons:

* FFT restrictions of flat-sky lensing codes often require that the cutout is square
* The cutouts returned will not actually have all side lengths of `boxLength` if we don't do this rotation, which the user explicitly requested

We want to express the positions of  all of our target lightcone objects in spherical coordinates, to perform the cutouts. Considering a single one of those objects; we want that coordinate system to be rotated such that the object of interest (which we will from now on assume is a halo) lies on the equator at

<p align="center">
(<i>r</i>, <i>&#x03B8;</i>, <i>&#x03D5;</i>) = (<i>r</i><sub>0</sub>, 90&#x00B0;, 0&#x00B0;)
</p>

where 

<p align="center">
<i>r</i><sub>0</sub> = (<i>x</i><sub><i>0</i></sub><sup>2</sup> + <i>y</i><sub><i>0</i></sub><sup>2</sup> + <i>z</i><sub><i>0</i></sub><sup>2</sup>)<sup>1/2</sup>.
</p>

Let's call the position vector of the halo before this rotation 

<p align="center">
<b>a</b> = [<i>x&#x2080;<i>, </i>y&#x2080;<i>, </i>z&#x2080;</i>], 
</p>

and after,

<p align="center">
 <b>b</b> = [<i>x</i><sub>rot</sub>, <i>y</i><sub>rot</sub>, <i>z</i><sub>rot</sub>] = [<i>r</i><sub>0</sub>, 0, 0]
</p>

We perform this rotation for the target halo via the *Rodrigues rotation formula*, which answers the following question: given a position vector **v**, a normalized axis of rotation **k**, and an angle of rotation *&#x03B2;*, what is an analytical form for a new vector **v**<sub>rot</sub> which is **v** rotated by an angle *&#x03B2;* about **k**?

First, we find **k** by taking the cross product of two vectors defining the 
plane of rotation. The obvious choice of these two vectors are **a** and **b**, as 
defined above;

<p align="center">
<b>k</b> = (<b>a</b> &#x2A2F; <b>b</b>) / &#x2016;<b>a</b> &#x2A2F; <b>b</b>&#x2016;
</p>

then, for any other position vector **v**, **v**<sub>rot</sub> is given by

<p align="center">
<b>v</b><sub>rot</sub> = <b>Rv</b>
</p>

where **R** is the rotation matrix through an angle *&#x03B2;* about the axis **k**, given as 

<p align="center">
<b>R</b> = <b>I</b> + sin(<i>&#x03B2;</i>)<b>K</b> + (1-cos(<i>&#x03B2;</i>))<b>K</b><sup>2</sup>.
</p>

**I**, here, is the 3x3 identity matrix, and **K** is the "cross-product matrix" for the unit vector **k**:

<p align="center">
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#x2308;&nbsp;&nbsp;&nbsp;0&nbsp;&nbsp;&nbsp;&nbsp;-<i>k</i><sub><i>z</i></sub>&nbsp;&nbsp;&nbsp;&nbsp;<i>k</i><sub><i>y</i></sub>&nbsp;&nbsp;&#x2309;<br>&nbsp;&nbsp;
<b>K</b>&nbsp;=&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;<i>k</i><sub><i>z</i></sub>&nbsp;&nbsp;&nbsp;&nbsp;0&nbsp;&nbsp;&nbsp;&nbsp;-<i>k</i><sub><i>x</i></sub>&nbsp;&nbsp;&nbsp;|&nbsp;<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#x230A;&nbsp;&nbsp;-<i>k</i><sub><i>y</i></sub>&nbsp;&nbsp;&nbsp;&nbsp;<i>k</i><sub><i>z</i></sub>&nbsp;&nbsp;&nbsp;&nbsp;0&nbsp;&nbsp;&nbsp;&#x230b;<br>
</p>

So for each target halo, **k**, *&#x03B2;*, and **R** are each computed once, and **v**<sub>rot</sub> is computed for each object (e.g. simulation particle), from the input lightcone, that should fill the resultant cutout. If we'd like to cut out *M* target halos from a particular lightcone, and have *N* particles filling that lightcone, then that means computing **v**<sub>rot</sub> *M*\**N* times, which may be a very large number. For efficiency, then, we follow this procedure (gorey details found in code comments):

1. Define the geometry of the cutout at the equator of our spherical coordinate system/on the *x*-axis of the base simulation's cartesian coordinate system. This means that *d&#x03B8;* = *d&#x03D5;* = *B*/2, where *B* is the `--boxLength` argument.

2. Define the vectors pointing to the four corners of the square bounded by *d&#x03B8;* and *d&#x03D5;* as **A**, **B**, **C**, and **D**.

3. Move the square field of view to the position of the target halo by inverting the rotation matrix:<br>**A**<sub>rot</sub> = **R**<sup>-1</sup>**A**<br>and similarly for **B**, **C**, and **D**. 

4. Make an "initial guess" around the field of view by doing a cut in constant *&#x03B8;* and *&#x03D5;* bounds, given by the maximum and minimum angular coordinates of **A**, **B**, **C**, and **D**. Add a 10 arcmin buffer.

5. For all lightcone objects (e.g. particles) surviving this initial cut, perform the proper rotation **v**<sub>rot</sub> = **Rv**. This number of objects will surely be &#x226A;*N*

If the 5 steps above don't make much sense, please let me know (contact info below) and perhaps I can put together an explanatory animation.

</p>
</details>

## Additional Arguments and Options

`-m` or `--massDef` specified the mass definition, in the case that `-f` is being used. Can either be `sod` or `fof`. If `sod`, then it is expected that each row of the file pointed to by `-f` has 7 elements, with teh fourth being the spherical overdensity radius. if `fof`, then the radius should be omitted from that file, else everything will break :)

`-v` or `--verbose` tells the application to generate tons of output, including explicity printing the rotation matrices and similar objects being used for Use Case 2

`--timeit` instruct the application to report wall-times for the data read, redistribution, cutout computation, and write-out

`--overwrite` allows the program to delete any contents inside of the `output directory`, rather than crashing with a warning

`--posOnly` will cause only particle positions, ids, and redshifts to be output. Velocities and lightcone rotation/replication information will be discarded (this should speed up both the data redistribution and the write-out).

For example, to run multiple cutouts under Use Case 2 with an `fof`mass definition, and turn on all of the options above, one would execute

```
lc_cutout <input lightcone directory> <output directory> <min redshift> <max redshift> --haloFile <input object file> --boxLength <box length> --massDef fof --verbose --tiemit --overwrite --posOnly
```

## Caveats

* Note that the Use Case 1 does not perform the coordinate rotation which is described in Use Case 2 (under the "click here to expand" details). So, cutouts returned will not necessarily be square, or symmetrical, if far from the coordinate equator. Even given Use Case 2, cutouts will not necessarily be square (though they should always be symmetrical) if the opening angle of the cutout breaks the small-angle approximation.

* The parallelism in this application occurs *spatially*, not temporally. That is, the lightcone *volume* is decomposed across MPI ranks, which prallelizes the read in, computation, and write-out. But there is no parallelism in *redshift*-space, meaning that each lightcone "step" (portion of the lightcone volume originating from a particular simulation snapshot) are treated in serial. Further, if option `-f` is used as described under Use Case 2, then those multiple requested cutouts are also treated serially. 

* The requested `min redshift` and `max redshift` are converted to a simulation step number assuming a simulation run that included 500 total time steps, and began at a redshift of 200. At the moment, there is no way for the user to easily change this, other than modifying the calls to `getLCSteps()` in `src/main.cpp` and rebuilding (the default values and parameter names controlling this info can be seen in the `getLCSteps()` function declaration in `src/util.h`).


#  Authors

Joe Hollowed

Argonne National Laboratory, Cosmological Physics and Advanced Computing Group (CPAC), 2018

jphollowed@gmail.com

## Acknowledgments

This code is based on earlier template code written by CPAC's Steve Rangel. I also received very helpful development contribution from Patricia Larsen at CPAC. Thanks friends :)

