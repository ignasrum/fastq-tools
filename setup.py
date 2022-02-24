from setuptools import setup

from fastq_tools.version import __version__


setup(
    name='fastq-tools',
    version=__version__,
    description='',
    url='https://github.com/ignasrum/fastq-tools',
    author='Ignas Rumbavicius',
    license='MIT',
    packages=['fastq_tools/fastq_compare',
              'fastq_tools/fastq_mixshuffle',
              'fastq_tools/fastq_rename'],
    python_requires='>=3.7',
    install_requires=['argparse',
                      'pytest',
                      'pysam'],
    entry_points = {
        'console_scripts': ['fastq-compare=fastq_tools.fastq_compare.compare:main',
                            'fastq-mixshuffle=fastq_tools.fastq_mixshuffle.mixshuffle:main',
                            'fastq-rename=fastq_tools.fastq_rename.rename:main'],
    }
)
