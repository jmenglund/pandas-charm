from setuptools import setup, find_packages
from os.path import join, dirname


setup(
    name='pandas-charm',
    version=__import__('pandascharm').__version__,
    description=(
        'A small Python library for getting character matrices '
        '(alignments) into and out of pandas'),
    long_description=open(
        join(dirname(__file__), 'README.rst'), encoding='utf-8').read(),
    packages=find_packages(exclude=['docs', 'tests*']),
    py_modules=['pandascharm'],
    install_requires=['pandas>=0.16', 'numpy', 'dendropy>=4'],
    extras_require={'test': ['coverage', 'pytest', 'pytest-cov']},
    author='Markus Englund',
    author_email='jan.markus.englund@gmail.com',
    url='https://github.com/jmenglund/pandas-charm',
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3'],
    keywords=['alignment', 'biopython', 'DendroPy', 'pandas'],
)
