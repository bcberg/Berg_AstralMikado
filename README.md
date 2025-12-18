# Berg_AstralMikado
Contains code used for Berg & Allard, *Astral architecture can enhance mechanical strength of cytoskeletal networks by modulating percolation thresholds*, publication pending.

## Contents
There are two components to the repository.
Code pertaining to mechanical simulations in Cytosim appears in ``mechanics/``, and code used to study geometric properties of astral networks (e.g., dangling ends, percolation) appears in ``geometry/``.

### ``mechanics/`` directory
Code is organized by type:
* ``analysis/``: Matlab scripts and functions used to read and analyze the outputs of Cytosim simulations.
* ``displays/``: files that control the appearance of simulation snapshots from Cytosim.
* ``drivers/``: Bash scripts used to initiate and organize simulation runs. 
* ``runs/``: Cytosim ``config.cym`` files used to run individual simulations, organized by simulation routine.
* ``templates/``: Templates for Cytosim configuration files (extension ``.cym.tpl``); each template corresponds to a simulation routine (subdirectories of ``runs/``).

### ``geometry/`` directory
Contains functions and scripts used to compute percolation probabilites, analyze dangling ends, and count mechanically productive segments.

## Running a simulation
First, you must clone this repository to your machine. 
The bash drivers assume the location ``~/Documents/Berg_AstralMikado``, but this can be modified.

### ``mechanics/``
Mechanical simulations require downloading and compiling a customized version of Cytosim, available at https://github.com/bcberg/cytosim-bcb.
Make note of where you install Cytosim.
The bash drivers assume it is located at ``~/Documents/cytosim-bcb``, but this can be modified.

To run and analyze a simulation, follow these steps:

1. Download both repositories (``Berg_AstralMikado`` and ``cytosim-bcb``).
2. Compile Cytosim on your machine.
Instructions are available in ``cytosim-bcb/README.md``.
3. Select a simulation routine from ``mechanics/templates/``, named with format ``simulation_name.cym.tpl``.
Each template file indicates its purpose and specifies sweep parameters (in double square brackets ``[[ ]]``).
4. Open the appropriate bash driver script (in ``mechanics/drivers/``).
For most simulations, this is ``simul_stand.sh`` (simulations handled by this script are listed in a comment at the top of the file).
If there is a driver named ``run_simulation_name.sh``, use this instead.

    1. Verify that the locations of ``Berg_AstralMikado`` and ``cytosim-bcb`` are correctly specified at the top of the file.
    2. Adjust the ``nproc`` parameter as desired (e.g., ``nproc=8`` enables simulations to run on up to 8 CPU cores in parallel).

5. Execute the driver script in your terminal.
This will typically run for several hours.
6. Once all simulations have finished, check that every individual run folder (``mechanics/runs/run????``) has a ``initial_pos.txt`` and ``final_pos.txt`` output file.
(This amounts to checking that each individual run folder contains the same number of files.)
    * In the event Cytosim failed to generate files for particular runs, try modifying the random seed near the top of each individual ``config.cym`` file.
    Note that this must be done for **every** run corresponding to a particular network's force-displacement curve.
        * Example: If the template specified a sweep over 6 force values, then ``run0000`` through ``run0005`` contain simulations of one particular network's response to each individual force value.
        If ``run0002`` failed to generate either ``.txt`` file, you must change the random seed in each configuration file (``run0000/config.cym`` through ``run0005/config.cym``) to a new number, e.g. by adding 1 to the failed random seed.
        * It is important that you change the random seed to the same new number in each of these files so that the simulations still compute a force-displacement curve for a single network.
        * After replacing failed random seeds with new ones, you must re-run ``cytosim-bcb/bin/sim`` in each directory that contains a updated ``config.cym`` file.

7. Analyze the output files by running the corresponding Matlab script in ``mechanics/analysis/``.
Adjust file paths/directories in these scripts as needed.

### ``geometry/``
Geometric computations require an installation of MATLAB (The MathWorks, Inc.).
They may be run entirely within MATLAB or with MATLAB_mex code (requires compiling the mex code on your machine).
Adjust this by setting the ``useMEX`` argument to ``true`` (allows use of mex code) or ``false`` (uses only MATLAB code).
It is highly recommended to have MATLAB's Parallel Computing Toolbox installed to accelerate percolation computations.

The core function of this section of code is ``generateAstralNetwork.m``.
It is called by higher-level functions and scripts:

* ``getPercCurve`` and ``getAllPercCurves`` are functions that compute percolation probability data for, respectively, a particular astral number or a list of astral numbers.
    * For details on the types of percolation the code checks for, see ``percCheck.m``.
* ``danglingEnds.m`` (script) quantifies dangling end properties in astral networks.
* ``springsPerNode.m`` (script) counts the number of productive filament segments per node in astral networks

Equal angle aster percolation is handled by ``generateAstralNetwork_eqSpaced.m`` and ``getAllPercCurves_eqSpaced.m``.

Estimation of critical percolation densities is done in ``percCrit.m`` and ``percCrit_eqSpaced.m``.

Functions are documented thoroughly in their defining ``.m`` file.