from setuptools import setup, find_packages

setup(
    name="fisix-lib",
    version="0.1.0",
    description="Python utility libraries for physics calculations and simulations",
    author="lspini6-code",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "iminuit",
    ],
    python_requires=">=3.7",
    url="https://github.com/lspini6-code/fisix-lib",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)