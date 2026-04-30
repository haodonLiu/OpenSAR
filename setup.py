"""Setup configuration for CSAR package."""

from setuptools import setup, find_packages

setup(
    name="csar",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.10",
    install_requires=[
        "rdkit>=2024.3.0",
        "pandas>=2.0.0",
        "numpy>=1.24.0",
        "networkx>=3.0",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "ruff>=0.1.0",
            "mypy>=1.0.0",
            "types-openpyxl>=3.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "csar=src.main:main",
        ],
    },
)
