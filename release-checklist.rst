Release checklist
=================

Things to remember when making a new release of pandas-charm.

#.  Changes should be made to some branch other than master (a pull request should then be created before making the release).

#.  Update the release (version) numbers in *setup.py* and *pandascharm.py*.

#.  Make desirable changes to the code.

#.  Run tests with coverage report and PEP8 check:

    .. code-block::

        $ py.test -v --cov-report term-missing --cov predsim.py --pep8

#.  Update the documentation in *README.rst*.

#.  Update *CHANGELOG.rst*.

#.  Create pull request(s) with changes for the new release.

#.  Create the new release in GitHub.

#.  Create distributions and upload the files to `PyPI <https://pypi.python.org/pypi>`_.

    .. code-block::

        $ python setup.py bdist_wheel --universal
        $ python setup.py sdist
 