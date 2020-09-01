from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='OceanLab',
    version='0.0.1',
    description='Python functions for Physical Oceanography',
    long_description=readme(),
    url='https://github.com/iuryt/OceanLab',
    author='Iury T. Simoes-Sousa',
    author_email='simoesiury@gmail.com',
    license='MIT',
    py_modules=['OA','EOF','DYN'],
    package_dir={'': 'src'},
)
