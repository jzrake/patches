# `patches`


## A C++14 micro-library for AMR in computational gasdynamics


## Goal
This code provides a `Database` structure, which maintains a sparse collection of logically cartesian grid patches at various levels of refinement. Fetching guard zones from neighboring patches at different refinement levels is done automatically.


## Dependencies
Depends on the [`nd::array`](https://github.com/jzrake/ndarray) multi-dimensional array library. Single-header available [here](https://github.com/jzrake/ndarray/blob/master/include/ndarray.hpp).


## Basic usage
A trivial example in Python-like pseudo-code is shown below. Note the library is in C++, but it's easier to convey the concepts if we pretend it's Python for a moment.

```python
# Create a database:
database = Database(100, 100)

# Populate it with 2x2 patches at the base level:
database.insert(i=0, j=0, level=0, make_array(100, 100])
database.insert(i=0, j=1, level=0, make_array(100, 100])
database.insert(i=1, j=0, level=0, make_array(100, 100])
database.insert(i=1, j=1, level=0, make_array(100, 100])

# Erase a block:
database.erase(i=1, j=1, level=0)

# Replace it with 4 patches at level 1:
database.insert(i=2, j=2, level=1, make_array(100, 100))
database.insert(i=2, j=3, level=1, make_array(100, 100))
database.insert(i=3, j=2, level=1, make_array(100, 100))
database.insert(i=3, j=3, level=1, make_array(100, 100))

# Get data at index (2, 2) on level 1, with 2 guard zones copied, prolonged,
# or restricted as needed:
my_array = database.fetch(i=2, j=2, level=1, guard=2)

# Solve some equations on the grid patch (reducing its size by 2),
# and merge in the result with a given Runge-Kutta parameter:
database.commit(advance(my_array, dt), rk_parameter=0.5)
```


## Status
Currently, the library supports a few simple 2D hydro codes that require static mesh refinement. There is no 3D support yet (although adding it is trivial). User-defined prolongation and restriction operators will be added soon. Currently data can be stored within cells, at mesh vertices, and on cell faces. Support for storage on edges will be added soon. So-called flux registers can be emulated by storing data on faces (or edges if solving MHD equations with constrained transport).

Support for MPI applications is planned. This will work by adding a `Database::fetch_async` method which returns a `std::promise<Array>` instead of an `Array` directly. However the actual code to place and fulfill remote queries is outside the scope of this module, and will be delegated to a user callback.
