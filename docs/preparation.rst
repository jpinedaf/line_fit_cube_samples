###################
Preparation of Data
###################

Before starting the line fit we need to prepare the data. As an example, if 
different lines need to be combined later in the analysis, then my suggestion 
is to first prepare all the cubes to a common beam and grid.

Convolving to a common beam
============================
The following code shows how to convolve a cube to a common beam. It is 
implemented to work witht a list of files to be convolved. The code will 
convolve the cubes to the beam of the first cube in the list.

.. literalinclude:: ../samples/common_beam.py
   :linenos:


File regridding
===============

The regridding of cubes to a matched grid can be done by leveraging the 
options already available in the `spectral-cube <https://spectral-cube.readthedocs.io/en/latest/>`_
package. The following code shows how to regrid a cube to a common grid.

.. literalinclude:: ../samples/regrid_cubes.py
   :linenos:

Unit conversion and primary beam correction
===========================================

The convertion of interferometric data into Kelvin units and including 
the primary beam conversion can be performed with the following code.

.. literalinclude:: ../samples/convert_flux_pb.py
   :linenos:

