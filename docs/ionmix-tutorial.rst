Ionmix Tutorial
````````````````
Welcome to the tutorial.

Getting Started
~~~~~~~~~~~~~~~~~~~~
Make sure you should see a folder called::

    ionmix_dp/

Within should be::

    doc/  examples/  scripts/ src/ test/

Let's begin with an example. Ionmix should come ready for compilation, so
go into src/ and give a make command to generate ionmix_dp. Ionmix is now
compiled. Enter ionmix_dp/examples/. You should see a folder, polystyrene/,
containing polystyrene-imx-007.inp. This is the input data for ionmix.
Ionmix only will take information from a file called ionmxinp, so make
a link to the .inp file with::

ln -s polystyrene-imx-007.inp ionmxinp

followed by a ../../src/ionmix_dp to run ionmix.

The files generated should be cnrdeos and ionmxout. ionmxout contains human
readable information about the material. cnrdeos is currently the file flash
reads at runtime containing the same information as ionmxout but mostly
unintelligible to a human.

Customized Materials
~~~~~~~~~~~~~~~~~~~~~~~
For organizational purposes, you will probably want to make a folder in 
ionmix_db/ called materials/.  Within this folder, keep folders of the materials
you are interested, like polystyrene/ in examples/.

Say we are interested in a material called new_material. In
materials/new_material/ generate or copy an existing .inp file. It is
encouraged notion to make the file called::

new_material-imx-00n.inp

where n is a number I will figure out later. For a reference, look at
polystyrene-imx-007.inp.

A complete list and explanation of the parameters (group boundaries, control
switches, constants, input variables, etc..) can be found in
docs/MacFarlane-1987.pdf.

Note that isw(6) can be set to either 1 or 3. 1 means the plasma is in local
thermal equilibrium and should be used for high densities. 3 means it is not
in LTE and should be chosen for low densities.

