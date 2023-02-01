from setuptools import setup


setup(
    name='biosampleMeta',
    python_requires='>3.10.0',
    version='0.0.1',
    url='https://github.com/veg/biosampleMeta',
    download_url="https://github.com/veg/biosampleMeta/archive/v0.0.1.tar.gz",
    description='NCBI metadata scraping',
    author='Stephen D. Shank',
    author_email='sshank314@gmail.com',
    maintainer='Stephen D. Shank',
    maintainer_email='sshank314@gmail.com',
    install_requires=[
        'biopython>=1.76',
        'beautifulsoup4>=4',
        'xmltodict>=0.13.0'
    ],
    packages=['biosampleMeta'],
    entry_points={
        'console_scripts': [
            'efetch = biosampleMeta.cli.efetch:efetch_cli',
            'sra_scrape = biosampleMeta.cli.efetch:efetch_cli'
        ]
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
    ]
)
