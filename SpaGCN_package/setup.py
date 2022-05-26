import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SpaGCN", 
    version="1.2.5",
    author="Jian Hu",
    author_email="jianhu@pennmedicine.upenn.edu",
    description="SpaGCN: Integrating gene expression and histology to identify spatial domains and spatially variable genes using graph convolutional networks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jianhuupenn/SpaGCN",
    packages=setuptools.find_packages(),
    install_requires=["python-igraph","torch","pandas","numpy","scipy","scanpy","anndata","louvain","sklearn", "numba"],
    #install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
