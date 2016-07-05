pandas-charm
============

.. image:: https://travis-ci.org/jmenglund/pandas-charm.svg?branch=master
    :target: https://travis-ci.org/jmenglund/pandas-charm

.. image:: https://codecov.io/gh/jmenglund/pandas-charm/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/jmenglund/pandas-charm

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: https://opensource.org/licenses/MIT

.. image:: https://zenodo.org/badge/23107/jmenglund/pandas-charm.svg
    :target: https://zenodo.org/badge/latestdoi/23107/jmenglund/pandas-charm

|

pandas-charm is a small Python package (or library) for getting character 
matrices (alignments) into and out of `pandas <http://pandas.pydata.org>`_.
The intention of the package is to make pandas interoperable with 
other scientific packages that can be used for working with character 
matrices, like for example `BioPython <http://biopython.org>`_ and 
`Dendropy <http://dendropy.org>`_.

With pandas-charm, it is currently possible to convert between the 
following objects:

* BioPython MultipleSeqAlignment <-> pandas DataFrame
* DendroPy CharacterMatrix <-> pandas DataFrame

Source repository: `<https://github.com/jmenglund/pandas-charm>`_


The name
--------

pandas-charm got its name from the pandas library plus an acronym for
CHARacter Matrix.


Installation
------------

For most users, the easiest way is probably to install the latest version 
hosted on `PyPI <https://pypi.python.org/>`_:

.. code-block::

    $ pip install pandas-charm

The project is hosted at https://github.com/jmenglund/pandas-charm and 
can be installed using git:

.. code-block::

    $ git clone https://github.com/jmenglund/pandas-charm.git
    $ cd pandas-charm
    $ python setup.py install


Running tests
-------------

After installing the pandas-charm, you may want to check that everything
works as expected. Below is an example of how to run the tests with pytest. 
The packages BioPython, DendroPy, pytest, coverage, and pytest-cov need 
to be installed.

.. code-block::

    $ cd pandas-charm
    $ py.test -v --cov-report term-missing --cov pandascharm.py


Usage
-----

Below are a few examples on how to use pandas-charm. The examples are 
written with Python 3 code, but pandas-charm should work also with 
Python 2.7. You need to install BioPython and/or DendroPy manually 
before you start:

.. code-block::

    $ pip install biopython
    $ pip install dendropy


Converting a DendroPy CharacterMatrix to a pandas DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: pycon

    >>> import pandas as pd
    >>> import pandascharm as pc
    >>> import dendropy
    >>> dna_string = '3 5\nt1  TCCAA\nt2  TGCAA\nt3  TG-AA\n'
    >>> print(dna_string)
    3 5
    t1  TCCAA
    t2  TGCAA
    t3  TG-AA
    
    >>> matrix = dendropy.DnaCharacterMatrix.get_from_string(
    ...     dna_string, schema='phylip')
    >>> df = pc.from_charmatrix(matrix)
    >>> df
      t1 t2 t3
    0  T  T  T
    1  C  G  G
    2  C  C  -
    3  A  A  A
    4  A  A  A

As seen above, characters are stored as rows and sequences as 
columns in the DataFrame. If you want rows to hold sequences, 
it is easy to transpose the matrix in pandas:

.. code-block:: pycon

    >>> df.transpose()
        0  1  2  3  4
    t1  T  C  C  A  A
    t2  T  G  C  A  A
    t3  T  G  -  A  A


Converting a pandas DataFrame to a Dendropy CharacterMatrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: pycon

    >>> import pandas as pd
    >>> import pandascharm as pc
    >>> import dendropy
    >>> df = pd.DataFrame({
    ...     't1': ['T', 'C', 'C', 'A', 'A'],
    ...     't2': ['T', 'G', 'C', 'A', 'A'],
    ...     't3': ['T', 'G', '-', 'A', 'A']})
    >>> df
      t1 t2 t3
    0  T  T  T
    1  C  G  G
    2  C  C  -
    3  A  A  A
    4  A  A  A
    
    >>> matrix = pc.to_charmatrix(df, type='dna')
    >>> print(matrix.as_string('phylip'))
    3 5
    t1  TCCAA
    t2  TGCAA
    t3  TG-AA


Converting a BioPython MultipleSeqAlignment to a pandas DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: pycon

    >>> from io import StringIO
    >>> import pandas as pd
    >>> import pandascharm as pc
    >>> from Bio import AlignIO
    >>> dna_string = '3 5\nt1  TCCAA\nt2  TGCAA\nt3  TG-AA\n'
    >>> f = StringIO(dna_string)  # make the string a file-like object
    >>> alignment = AlignIO.read(f, 'phylip-relaxed')
    >>> print(alignment)
    SingleLetterAlphabet() alignment with 3 rows and 5 columns
    TCCAA t1
    TGCAA t2
    TG-AA t3
    >>> df = pc.from_bioalignment(alignment)
    >>> df
      t1 t2 t3
    0  T  T  T
    1  C  G  G
    2  C  C  -
    3  A  A  A
    4  A  A  A


Converting a pandas DataFrame to a BioPython MultipleSeqAlignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: pycon

    >>> import pandas as pd
    >>> import pandascharm as pc
    >>> import Bio
    >>> df = pd.DataFrame({
    ...     't1': ['T', 'C', 'C', 'A', 'A'],
    ...     't2': ['T', 'G', 'C', 'A', 'A'],
    ...     't3': ['T', 'G', '-', 'A', 'A']})
    >>> df
      t1 t2 t3
    0  T  T  T
    1  C  G  G
    2  C  C  -
    3  A  A  A
    4  A  A  A
    
    >>> alignment = pc.to_bioalignment(df, alphabet='generic_dna')
    >>> print(alignment)
    SingleLetterAlphabet() alignment with 3 rows and 5 columns
    TCCAA t1
    TGCAA t2
    TG-AA t3
    

License
-------

pandas-charm is distributed under 
`the MIT license <https://opensource.org/licenses/MIT>`_.


Citing
------

If you use results produced with this package in a scientific 
publication, please just mention the package name in the text and 
cite the Zenodo DOI of this project:

.. image:: https://zenodo.org/badge/23107/jmenglund/pandas-charm.svg
    :target: https://zenodo.org/badge/latestdoi/23107/jmenglund/pandas-charm

You can select a citation style from the dropdown menu in the 
*"Cite as"* section on the Zenodo page.
