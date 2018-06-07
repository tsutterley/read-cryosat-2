from setuptools import setup, find_packages
setup(
	name='read-cryosat-2',
	version='1.0.0.3',
	description='Reads and writes data from the ESA CryoSat-2 mission',
	url='https://github.com/tsutterley/read-cryosat-2',
	author='Tyler Sutterley',
	author_email='tyler.c.sutterley@nasa.gov',
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
	install_requires=['numpy','h5py','future'],
)
