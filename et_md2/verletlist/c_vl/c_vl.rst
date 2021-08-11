This file documents a binary extension module built from C++ code.

You should document the Python interfaces, *NOT* the C++ interfaces.

Module et_md2.verletlist.c_vl
*****************************
This is a C++ binary extension. It links a shared library vl_lib.

Module :py:mod:`c_vl` built from C++ code in :file:`et_md2/verletlist/c_vl/c_vl.cpp`.

.. function:: add(x,y,z)
   :module: et_md2.verletlist.c_vl
   
   Compute the sum of *x* and *y* and store the result in *z* (overwrite).

   :param x: 1D Numpy array with ``dtype=numpy.float64`` (input)
   :param y: 1D Numpy array with ``dtype=numpy.float64`` (input)
   :param z: 1D Numpy array with ``dtype=numpy.float64`` (output)
   