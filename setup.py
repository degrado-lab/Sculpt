from setuptools import setup, find_packages
from pathlib import Path

# Read the long description from README.md
long_description = Path(__file__).with_name("README.md").read_text(encoding="utf8")

setup(
    name="sculpt",
    version="0.1.0",
    description="Sculpt — geometry-first, iterative enzyme design",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Nicholas Freitas / degrado-lab",
    url="https://github.com/degrado-lab/Sculpt",
    packages=find_packages(exclude=("tests", "__pycache__")),
    include_package_data=True,
    install_requires=[
        # Ribbon is required by the library (see README). Pin a permissive minimum.
        "ribbon-toolkit>=0.3.1",
        "mdtraj",
    ],
    python_requires=">=3.12",
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
