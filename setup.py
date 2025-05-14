from setuptools import setup, find_packages


setup(
    name='CRISPRlungo',
    version='0.1',
    package_dir={'': 'src'},             
	py_modules=['CRISPRlungo', 
		'CRISPRlungo_custom_mutation', 
		'CRISPRlungo_umi', 
		'CRISPRlungo_insert_analysis', 
		'CRISPRlungo_regular', 
		'CRISPRlungo_visualization', 
		'CRISPRlungo_minimap', 
		'CRISPRlungo_single_map'],   
    entry_points={
        'console_scripts': [
            'CRISPRlungo = CRISPRlungo:main'  
        ],
    },
)
