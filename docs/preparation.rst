###################
Preparation of Data
###################

Before starting the line fit we need to prepare the data. As an example, if 
different lines need to be combined later in the analysis, then my suggestion 
is to first prepare all the cubes to a common beam and grid.


File regridding
===============

I will add a sample file showing how to modify a cube to match 
angular resolution and pixel grid to that of a second line.


Unit conversion and primary beam correction
===========================================

The convertion of interferometric data into Kelvin units and including 
the primary beam conversion can be performed with the following code.

.. literalinclude:: ../samples/convert_flux_pb.py
   :linenos:

