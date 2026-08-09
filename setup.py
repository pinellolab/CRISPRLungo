from setuptools import setup, find_packages


setup(
    name='CRISPRlungo',
    version='0.3',
    package_dir={'': 'src'},      
    packages=find_packages(where='src'),
    include_package_data=True,
	py_modules=['CRISPRlungo', 
		'CRISPRlungo_log',
		'CRISPRlungo_umi', 
		'CRISPRlungo_insert_analysis', 
		'CRISPRlungo_visualization', 
		'CRISPRlungo_minimap', 
		'CRISPRlungo_single_map',
        'CRISPRlungo_mutation_analysis',
        'CRISPRlungoAllele',
        'CRISPRlungoBatch'],   
    package_data={
        "CRISPRlungo_assets": [
            "css/*"
            "css/bootstrap-5.3.7-dist/*",
            "css/bootstrap-5.3.7-dist/css/*",
            "css/bootstrap-5.3.7-dist/js/*",
            "css/bootstrap-5.3.7-dist/fonts/*",
            "css/bootstrap-5.3.7-dist/**/*",
        ]
    },
    entry_points={
        'console_scripts': [
            'CRISPRlungo = CRISPRlungo:main',
            'CRISPRlungoAllele = CRISPRlungoAllele:main',
            'CRISPRlungoBatch = CRISPRlungoBatch:main'
        ],
    },
    
)
