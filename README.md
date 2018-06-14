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

Theta - coaltitutde coordinate
Phi - azimuthal coordinate
*x*, *y*, *z* - comoving cartesian coordinates
`lc_cutout` - the executable resultant from following the installation instructions

Two use cases are supported: 

### Use case 1: Constant angular bounds

Define the theta and phi bounds explicitly:
  
```
lc_cutout <input lightcone directory> <output directory> <depth> --theta <Theta_center> <dTheta> --phi <Phi_center> <dPhi>
```

where the Phi\_center argument is the azimuthal coordinate of the center of the field of view that one wishes to cut out of the lightcone, and dPhi is the angualar distance from this center to the edge of the cutout, and likewise for the similar Theta  args. That is, the result will be a sky area that spans 
 
```
(Theta_center - dTheta) < Theta < (Theta_center + dTheta)
(Phi_center - dPhi) < Phi < (Phi_center + Phi)
```

The expected angular units are DEGREES. The `--theta` and `--phi` flags can be replaced with `-t` and `-p`.

The `depth` parameter is the maximum desired redshift of the cutout (limited, of course, by the maximu redshift of the input lightcone catalog).


### Use case 2: Nonlinear angular bounds

Allow the Theta and Phi bounds to be computed interanally to obtain a cutout of a certain width (*box length*), in Mpc/h, centered on a certain cartesian positon, (*x*&#x2080;, *y*&#x2080;, *z*&#x2080;) Mpc/h (intended to be used for making cutouts centerd on specific simulation objects, like halos):

```
lc_cutout <input lightcone directory> <output directory> <depth> --halo <x_0> <y_0> <z_0> --boxLength <box length>
```

where the `--halo` and `--boxLength` flags can be replaced with `-h` and `-b`. Those options give the position of the object to center the cutout on, and the one-dimensional size of the field of view in Mpc/h, tangent to the line of sight, at the distance of the object of interest, respectively.

We want to express the positions of  all of our lightcone objects in spherical coordinates, to perform the cutout, and we want that coordinate system to be rotated such that the object of intererst lies on the equator at

(*r*, 90&#x00B0;, 0&#x00B0;)
where *r* = (*x*^2 + *y*^2 + *z*^2)&#x00BD;

Let's call the position vector of the halo before this rotation 
**a** = [*x*&#x2080;, *y*&#x2080;, *z*&#x2080;], and after, **b** = [*x*&#x1D63;&#x2092;&#x209C;, *y*&#x1D63;&#x2092;&#x209C;, *z*&#x1D63;&#x2092;&#x209C;] = [*r*, 0, 0]

We perform this rotation for each lightcone object via the *Rodrigues rotation formula*, which answers the following question: given a position vector **v**, a normalized axis of rotation **k**, and an angle of rotation &#x03B2;, what is an analytical form for a new vector **v**&#x1D63;&#x2092;&#x209C; which is **v** rotated by an anlge &#x03B2; about **k**?

First, we find **k** by taking the cross product of two vectors defining the 
plane of rotation. The obvious choice of these two vectors are **a** and **b**, as 
defined above;

k = (**a** &#x2A2F; **b**) / &#x2016;**a** &#x2A2F; **b**&#x2016;

then, for any other position vector **v**, **v**&#x1D63;&#x2092;&#x209C; is given by

**v**&#x1D63;&#x2092;&#x209C; = **v**cos&#x03B2; + (**k** &#x2A2F; **v**)sin&#x03B2; + **k**(**k** &#x22C5; **v**)(1-cos&#x03B2;)

This coordinate rotation is required because the bounding functions which define the field of view of the observer, while constant in theta-phi space, are nonlinear in cartesian space. The distortion is maximized near the poles of the spherical coordinate system, and minmized at the equator. Areas defined by constant theta-phi bounds then appear trapezoidal to the observer when far from the equator. It is important that our cutout areas are maintained as square for at least two reasons:

* At the moment, FFT restrictions of flat-sky lensing code require that the cutout is square
* The cutouts returned will not actually have all side lengths of `boxLength` if we don't do this rotation, which the user explicitly requested

Note that the first use case describe does not perform the coordinate rotation which is described in the second. So, cutouts returned will not necessarily be square, or symmetrical.

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

