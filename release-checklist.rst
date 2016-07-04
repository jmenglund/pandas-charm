Release checklist
=================

Things to remember when making a new release of pandas-charm.

1.  Changes should be made to the branch named "develop" (a pull request 
    should then be created before making the release).

2.  Update the release (version) number in *setup.py* and *pandascharm.py*.

3.  Make disirable changes to the code.

4.  Update the documentation in *README.rst*.

5.  Update *CHANGELOG.rst*.

6.  Create pull request(s) with changes for the new release.

7.  Create the release in GitHub.

8.  Create distributions and upload the files to PyPI

    .. code-block:: console
    
        $ python setup.py bdist_wheel --universal
        $ python setup.py sdist
