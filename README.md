# DalgGuessing

[![arXiv](https://img.shields.io/badge/arXiv-2510.26869-red.svg)](https://arxiv.org/abs/2510.26869)

## Introduction

This GitHub page is devoted to the Maple package **DalgGuessing**.

DalgGuessing is a standalone package for the guessing software from the Maple package [NLDE](https://github.com/T3gu1a/D-algebraic-functions).
The package currently exports two commands for recovering D-algebraic functions from the first few coefficients of their power series. These are:

- ``DalgFunGuess`` :
- ``modDalgFunGuess``:


## Requirements

A recent version of [Maple](https://www.maplesoft.com/) (version $\geq$ 2018).

## Installation

The package was created as a standard Maple package. To use it, the Maple variable 
``libname`` should contain the path to the directory where the ``DalgGuessing.mla`` file is located. 

```
libname := "PATH_TO/DalgFunGuess.mla", libname:
```

Alternatively, this file should be in the same directory as the Maple worksheet being used.
After completing this setup, one loads the package using the command

```
with(DalgGuessing);
```
```math
[\texttt{DalgFunGuess},\texttt{modDalgFunGuess}]
```
## Author

- [Bertrand Teguia Tabuguia](https://bertrandteguia.com), [Email Me](mailto:yourname@example.com)
- licence: GNU General Public Licence v3.0.
