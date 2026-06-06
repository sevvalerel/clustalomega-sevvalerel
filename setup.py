from setuptools import setup, find_packages

setup(
    name='clustalomega-sevvalerel',
    version='1.0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'clustalomega=clustalomega.cli:main',
        ],
    },
    python_requires='>=3.8',
    description='Clustal Omega MSA algoritması',
    author='Sevval Erel',
)