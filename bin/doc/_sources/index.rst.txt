=====================
Project Documentation
=====================

License
=======
Copyright © 2021 lascauje

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

Philosophy
==========
This project implements the *Quadratic Inverse Problems in Phase retrieval*
from the scientific paper
`Convex-Concave Backtracking
for Inertial Bregman Proximal Gradient Algorithms in Non-Convex Optimization
<https://paperswithcode.com/paper/convex-concave-backtracking-for-inertial>`_

I have discovered a truly well-written explanation of this,
which this 80 column is too narrow to contain.

Furthermore, the goals were both to implement a scientific paper
and to make a comparison between two programming languages: clojure and python.

**This may not be the best clojure (respectively, python) code,
but it is a clojure (respectively, python) code.**

Code Documentation
==================
For the documentation please visit `clojure-doc <clj/bpg.core.html>`_,
and `python-doc <py/reference.html>`_.

In most function documentation there is: *See p.x (y.z)*.
It means "see the formula y.z page x of the scientific paper in pdf format"

The BPG acronym stands for Bregman Proximal Gradient.

Tools
=====
This section describes the tools used for the clojure (respectively, python) project.

============= ============================== ========================
tool          clojure                        python
============= ============================== ========================
code (M-x)    cider-mode                     python-mode
algebra       neanderthal                    numpy
build         leiningen                      poetry, nox
lint          cljfmt, eastwood, kondo, kibit flake8, darglint, pytype
test          leiningen, clojure.test        pytest
doc           codox                          sphinx
============= ============================== ========================

How to Build It?
================
For the clojure project, the following command runs linters, tests
and produces the documentation.

.. code-block:: bash

   $ lein ci

Respectively, the following command is for the python project.

.. code-block:: bash

   $ nox

The following command builds this documentation.

.. code-block:: bash

   $ ./build_doc.sh

The following code block describes the system dependencies.
It can be used to build a reproducible environment.
For clojure dependencies (respectively, python), visit
the file *project.clj* (respectively, *pyproject.toml*).

.. code-block:: docker

   python=3.9.1
   pip=20.2.4
   poetry=1.1.5
   nox=2020.12.31
   leiningen=2.9.5
   openjdk=11.0.10

Benchmark
=========
All computations have been run on an (old) computer:
Intel(R) Core(TM) i7-4600U CPU @ 2.10GHz, 2 cores, 3GiB, GNU/Linux.

Matrices are of dimension 10, 32, 128, 256, 512, or 1024,
which are duplicated 100 times, and with 1000 iterations.
Matrices are generated through a homemade linear congruential generator.
This ensures that the calculations are the same between
the both clojure and python implementations.

The scientific paper uses only matrices of 10x10 dimension
with 1000 iterations.

The following table is a comparison about execution time.

Results are reported in milliseconds, the lower the better.

========= ========= =========
size      clojure   python
========= ========= =========
10x10     3579.37   4534.34
32x32     4353.38   4856.39
128x128   25305.41  19917.74
256x256   62457.68  55281.28
512x512   210053.55 201830.46
1024x1024 986261.71 966849.65
========= ========= =========

The following table is a comparison about memory usage.

Results are reported in MB, the lower the better.

NB: results below are expressed **without** the repl memory footprint.
The clojure repl consumes about 500MB, and the python one 22MB.

========= ========= =========
size      clojure   python
========= ========= =========
10x10     42        1
32x32     73        3
128x128   320       12
256x256   401       55
512x512   835       203
1024x1024 844       770
========= ========= =========

How to Benchmark It?
~~~~~~~~~~~~~~~~~~~~
For the clojure execution time benchmark,
the following commands have been used.
*N* represents the matrix dimension, here 10, 32, 128, 256, 512 and 1024.

.. code-block:: bash

   $ lein repl

.. code-block:: clojure

   bpg.core=> (bench N)

Respectively, for the python benchmark.

.. code-block:: bash

   $ poetry run python

.. code-block:: python

   >>> from bpg import bpg

.. code-block:: python

   >>> bpg.bench(N)

For the memory usage benchmark, the *htop* command line has been used
to monitor manually clojure (respectively, python) process
from the repl start to the end of the computation.

A Salute
========
The repl driven development in clojure was a little bit of honey
thanks to cider-mode, and I think that
the python ecosystem is more mature (tools, documentations, etc.)

The design was first implemented in clojure, and then in python.

The first benchmark is quite about neanderthal and numpy framework, while
the second is about clojure (respectively, python) core language.

Even though the execution time (respectively, memory usage)
is very similar between the both languages,
python remains the fastest and consumes less memory than clojure.

NB: neanderthal is *work in progress*, with 13 contributors,
whereas numpy is *stable version* with more than 1,100 contributors.

Furthermore, the clojure repl memory footprint is **heavy**!
It's partly because it starts a local server to allow
interactive development (repl driven development, live debugging, etc.)
But ok, it's always too much. Could graalvm do better?

Any suggestion is welcome, feel free to open a pull request.

References
==========
- M.C. Mukkamala and P. Ochs and T. Pock and S. Sabach. 2019.
  `Convex-Concave Backtracking
  for Inertial Bregman Proximal Gradient Algorithms in Non-Convex Optimization
  <https://paperswithcode.com/paper/convex-concave-backtracking-for-inertial>`_

- Daniel Higginbotham. 2015.
  `Clojure for the Brave and True <https://www.braveclojure.com/>`_

- Dragan Djuric. 2021.
  `Neanderthal, Fast native-speed matrix and linear algebra in Clojure
  <https://neanderthal.uncomplicate.org/>`_

- Claudio Jolowicz. 2020.
  `Hypermodern Python
  <https://cjolowicz.github.io/posts/hypermodern-python-01-setup/>`_

- al. 2021.
  `Numpy, The fundamental package for scientific computing with Python
  <https://numpy.org/>`_
