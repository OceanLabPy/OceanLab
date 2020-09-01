from setuptools import setup


with open('README.md','r') as fh:
    long_description = fh.read()

setup(
    name='OceanLab',
    version='0.0.1',
    description='Python functions for Physical Oceanography',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/iuryt/OceanLab',
    author='Iury T. Simoes-Sousa',
    author_email='simoesiury@gmail.com',
    license='MIT',
    py_modules=['OA','EOF','DYN'],
    package_dir={'': 'src'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = [
        'seawater ~= 3.3',
        'numpy ~= 1.18',
    ],
)
