from setuptools import setup, find_packages


setup(
    name='CRISPRlungo',
    version='0.1',
    package_dir={'': 'src'},      
    packages=find_packages(where='src'),
    include_package_data=True,
	py_modules=['CRISPRlungo', 
		'CRISPRlungo_custom_mutation', 
		'CRISPRlungo_umi', 
		'CRISPRlungo_insert_analysis', 
		'CRISPRlungo_regular', 
		'CRISPRlungo_visualization', 
		'CRISPRlungo_minimap', 
		'CRISPRlungo_single_map',
        'CRISPRlungoAllele'],   
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
            'CRISPRlungoAllele = CRISPRlungoAllele:main'
        ],
    },
    
)
