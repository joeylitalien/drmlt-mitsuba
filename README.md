# Delayed Rejection Metropolis Light Transport

This repository contains the *bold-then-timid* implementation of [Delayed Rejection Metropolis Light Transport](https://joeylitalien.github.io/publications/drmlt) based on the [Mitsuba v0.6](https://www.mitsuba-renderer.org/download.html) renderer. Note that this work is a fork of the SIGGRAPH 2018 course [Light Transport Simulation in the Gradient Domain](https://github.com/beltegeuse/gradient-mts) repository.

If you want to understand the algorithm by looking at the code, you should start with:

  - Multiplexed MLT__
    - `./include/mitsuba/bidir/path_sampler.h`
    - `./src/libbidir/path_sampler.cpp`
  - __Delayed Rejection MLT__
    -  `./src/integrator/drmlt/*`

In case of problems/questions/comments don't hesitate to contact the authors directly.
 
Dependencies
------------

### Mitsuba
- [Boost](https://www.boost.org/)
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [OpenEXR](https://www.openexr.com/)
- [FFTW](http://www.fftw.org/)
- [LibPNG](http://www.libpng.org/pub/png/libpng.html)
- [ZLib](https://zlib.net/)
- [GLEW](http://glew.sourceforge.net/)
- [OpenImageIO](https://github.com/OpenImageIO/oiio)

### Python
- [PyEXR](https://github.com/tvogels/pyexr)
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)

Compiling
---------

We provided installation instructions for ArchLinux only, but the installation procedure is similar on Ubuntu using `apt install`. To install Mitsuba's dependencies, run:
```bash
pacman -Sy make boost eigen gcc openexr python3 fftw libpng jasper zlib cmake git awk xerces-c xorg glew openimageio python-pip
```

Configure and build the project using [Cmake](https://cmake.org/): 
```bash
mkdir build && cd "$_"
cmake ../
make -j
```

To install the tooling dependencies, run:
```bash
pacman -Sy python3 python-pip \
pip install numpy matplotlib pyexr
```

Delayed Rejection Framework
---------------------------

Our implementation of delayed rejection support three types of frameworks:
 - __Original Framework:__ Proposed by [Tierney & Mira [1999]](https://www.researchgate.net/publication/2767014_Some_Adaptive_Monte_Carlo_Methods_for_Bayesian_Inference). Suffers from vanishing acceptance at the second stage.
 - __Generalized Framework:__ Proposed by [Green & Mira [2001]](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.20.7698&rep=rep1&type=pdf). Uses reversible jump MCMC to solve the vanishing acceptance problematic by sampling an expensive intermediate state. 
 - __Pairwise Orbital:__ Based on the original framework but use an orbital mutations strategy at the second stage to solve the vanishing acceptance problem without the extra sampling computational overhead.


Path Sampling Technique
-----------------------

An important change from the previoux Mitsuba implementation is that now integrator can be run over three different path sampling techniques:
- __Unidirectional Path Tracing (PT):__ Unidirectional volumetric path tracer.
- __Bidirectional Path Tracing (BDPT):__ Bidirectional path tracer with Multiple Importance Sampling (MIS).
- __Multiplexed MLT (MMLT):__ [Multiplexed MLT](https://www.ci.i.u-tokyo.ac.jp/~hachisuka/mmlt.pdf) Bidirectional path tracer.

You can find these in `src/libbidir/path_sampler.cpp` and `include/bidir/path_sampler.h`.

Integrators
-----------

We modified and added the following integrators:
 - `src/integrator/pssmlt`: This is a modified version of the original PSSMLT algorithm.
 - `src/integrator/drmlt`: This is the core of our DRMLT algorithm.

You can select between them by using the `-D integrator=[pssmlt,drmlt]` command-line parameter.

### `pssmlt`
| Parameter | Description | Requirement |
|:----------|:------------|:--|
| `technique` | Path sampling technique | Required  (Options: `path, bdpt, mmlt`) |
| `kelemenStyleMutation` | Use Kelemen or Gaussian mutation | Optional (Default: `true`)  |
| `mutationSizeLow` | Kelemen lower bound | Optional (Default: `1/1024`) |
| `mutationSizeHigh` | Kelemen higher bound | Optional (Default: `1/64`) |
| `sigma` | Standard deviation of Gaussian mutation | Optional (Default: `1/64`) |

For example:
```bash
<PATH_TO_MITSUBA_BIN>/mitsuba <PATH_TO_SCENE>/scene.xml -D integrator=pssmlt -D technique=path                      
```

### `drmlt`

| Parameter | Description | Requirement |
|:----------|:------------|:--|
| `technique` | Path sampling technique | Required  (Options: `path, bdpt, mmlt`) |
| `type` | Delayed rejection framework | Required  (Options: `mira, green, orbital`) |
| `acceptanceMap` | Output acceptance map | Optional (Default: `false`)  |
| `timidAfterLarge` | Perform second stage after a large step | Optional (Default: `false`)  |
| `fixEmitterPath` | Fix emitter subpath during the second stage (Only with the `mmlt` technique) | Optional (Default: `false`)  |
| `useMixture` | Use an equal weight mixture of both stage and regular Metropolis-Hastings instead of DR | Optional (Default: `false`)  |
| `sigma` | Standard deviation of Gaussian mutation | Optional (Default: `1/64`) |
| `scaleSecond` | Scaling ratio of the second stage mutation | Optional (Default: `1/10`) |

For example:
```bash
<PATH_TO_MITSUBA_BIN>/mitsuba <PATH_TO_SCENE>/scene.xml \
                              -D integrator=drmlt       \
                              -D technique=mmlt         \
                              -D type=orbital           \
                              -D fixEmitterPath=true    \
                              -D acceptanceMap=false
```


Acceptance Map
--------------

When using the `drmlt` integrator, you can generate an acceptance map using the `-D acceptanceMap=true` option. Doing so will generate an RGB image such that the _R_-channel corresponds to the number of accepted samples at the first stage and the _G_-channel is the same for the second stage. To convert this image to a heatmap, use the standalone script `./tools/stages_heatmap.py`. For example, the following command saves the acceptance map during rendering:

```bash
<PATH_TO_MITSUBA_BIN>/mitsuba <PATH_TO_SCENE>/scene.xml \
                              -D integrator=drmlt       \
                              -D technique=bdpt         \
                              -D type=orbital           \
                              -D acceptanceMap=true
```

To generate the actual heatmap, run:

```bash
python <PATH_TO_MITSUBA_ROOT>/tools/stages_heatmap.py    \
          -t <PATH_TO_ACCEPTANCE_MAP>/acceptance_map.exr \ 
          -c [0.2,0.8]
```

| Parameter | Description | Requirement |
|:----------|:------------|:--|
| `t` | Acceptance map | Required |
| `c` | Pixel range (clip) for heatmap images | Optional (Default: `[0,1]`) |


Mixture
-------

To generate a comparison of our method against a naïve mixture of both stage, use the `-D useMixture=true` option under the `drmlt` integrator.


Scenes
------

- [Swimming Pool](http://beltegeuse.s3-website-ap-northeast-1.amazonaws.com/research/2020_DRMLT/scenes/swimming-pool_pssmlt.zip)
- [Aquarium](http://beltegeuse.s3-website-ap-northeast-1.amazonaws.com/research/2020_DRMLT/scenes/aquarium_mmlt.zip)
- [Veach Door](http://beltegeuse.s3-website-ap-northeast-1.amazonaws.com/research/2020_DRMLT/scenes/veach-door_mmlt.zip)
- [Glass of Water](http://beltegeuse.s3-website-ap-northeast-1.amazonaws.com/research/2020_DRMLT/scenes/glass-of-water_pssmlt.zip)


Change Logs
-----------

- 2020/07/29: Initial code release


License
-------

This code is released under the GNU General Public License (version 3).

This source code includes the following open source implementations:

- Screened Poisson reconstruction code from NVIDIA, released under the new BSD license.
- Mitsuba 0.6.0 by Wenzel Jakob, released under the GNU General Public License (version 3).
- A small part of [Tungsten](https://github.com/tunabrain/tungsten) by Benedikt Bitterli.
