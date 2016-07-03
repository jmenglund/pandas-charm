pandas-charm
============

pandas-charm is a small Python library for getting character matrices
(alignments) into and out of `pandas <http://pandas.pydata.org>`_.
Use the library to make pandas interoperable with
`BioPython <http://biopython.org>`_ and `Dendropy <http://dendropy.org>`_.

With pandas-charm, you can convert the following objects:

* BioPython MultipleSeqAlignment to/from pandas DataFrame
* DendroPy CharacterMatrix to/from pandas DataFrame


Source repository: `<https://github.com/jmenglund/pandas-charm>`_

.. image:: https://travis-ci.org/jmenglund/pandas-charm.svg?branch=master
    :target: https://travis-ci.org/jmenglund/pandas-charm

.. image:: https://codecov.io/gh/jmenglund/pandas-charm/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/jmenglund/pandas-charm


The name
--------

pandas-charm got its name from the library pandas plus an acronym for
CHARacter Matrix.


Installation
------------

The project is hosted at https://github.com/jmenglund/pandas-charm and 
can be installed using git:

.. code-block:: console

    $ git clone https://github.com/jmenglund/pandas-charm.git
    $ cd pandas-charm
    $ python setup.py install


Running tests
-------------

After installing the library, you may want to check that everething
works as expected. Below is an example of how to run the tests. The packages
BioPython, DendroPy, pytest, coverage, and pytest-cov need to be installed.

.. code-block:: console

    $ cd pandas-charm
    $ py.test -v --cov-report term-missing --cov pandascharm.py


Usage
-----

Below are a few examples (for Python 3) on how to use the library. 
Sequences are treated as rows and characters as columns in the DataFrame.
You need to install BioPython and DendroPy manually before you start:

.. code-block:: console

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
    
