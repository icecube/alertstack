import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="alertstack",
    version="0.0.1",
    author="Robert Stein",
    author_email="robert.stein@desy.de",
    description="Package for unbinned likelihood analysis of high-purity physics data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords="astroparticle physics science unbinned likelihood IceCube",
    url="https://github.com/icecube/alertstack",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "numpy",
        "healpy",
        "scipy",
        "matplotlib",
        "astropy",
        "coveralls"
    ],
    package_data={'alertstack': [
        'alertstack/data/*']},
    include_package_data=True
)

