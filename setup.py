from setuptools import setup, find_packages
setup(
	name='read-cryosat-2',
	version='1.0.0.8',
	description='Reads and writes data from the ESA CryoSat-2 mission',
	url='https://github.com/tsutterley/read-cryosat-2',
	author='Tyler Sutterley',
	author_email='tsutterl@uw.edu',
	license='MIT',
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Physics',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2',
		'Programming Language :: Python :: 2.7',
	],
	keywords='CryoSat-2 radar altimetry SIRAL',
	packages=find_packages(),
	install_requires=['numpy','h5py','netCDF4','future','lxml'],
)
