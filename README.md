
# Minimum k-Vertex Connected Graph Search

This repository implements the minimum k-VC search algorithm **VCtoB** proposed in our paper. If you find this code useful in your research, please cite our paper.

## Table of Contents
- [Minimum k-Vertex Connected Graph Search](#minimum-k-vertex-connected-graph-search)
  - [Table of Contents](#table-of-contents)
  - [Prerequisites](#prerequisites)
  - [Compilation](#compilation)
  - [Datasets](#datasets)
  - [Data Format and Conversion](#data-format-and-conversion)
  - [Usage](#usage)
    - [Parameters:](#parameters)
  - [Example](#example)

## Prerequisites
Ensure you have the following installed:
- C++ Compiler (supports C++11)

## Compilation
To compile the code, use the following command in the terminal:

```sh
g++ -std=c++11 -O3 -o VCtoB main.cpp Graph.cpp sbundle_tool.cpp
```

This command will produce an executable named `VCtoB`, which implements our algorithm.

## Datasets
The datasets used in our paper can be obtained from the following source:
- [139 Real-world Graphs Collection](http://lcs.ios.ac.cn/~caisw/Resource/realworld%20graphs.tar.gz)

## Data Format and Conversion
Our program accepts binary input data files only (`*.bin`). If your data files are in a different format, you can convert them using the `translate.cpp` tool included in the `data_transform_tool` directory. For example, to convert a `.mtx` file to the binary format (`.bin`), you can use the following commands:

```sh
g++ -g -o translate translate.cpp
./translate bio-celegans.mtx bio-celegans.bin
```

## Usage
To run the `VCtoB` executable, use the following syntax:

```sh
./VCtoB -g {path_to_graph} -k {k_value}
```

### Parameters:
- `-g` specifies the path to the graph file.
- `-k` specifies the k-value for the vertex connectivity.

## Example
Here is an example command to compute the exact maximum 6-bundle for the 'tech-as-caida2007.bin' graph with k=16:

```sh
./VCtoB -g datasets/tech-as-caida2007.bin -k 16
```