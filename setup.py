from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='read-cryosat-2',
    version='1.0.1.5',
    description='Reads and writes data from the ESA CryoSat-2 mission',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/tsutterley/read-cryosat-2',
    author='Tyler Sutterley',
    author_email='tsutterl@uw.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='CryoSat-2 radar altimetry SIRAL',
    packages=find_packages(),
    install_requires=['numpy','scipy','h5py','netCDF4','gdal','fiona',
        'geopandas','shapely','pyproj','future','lxml'],
)
