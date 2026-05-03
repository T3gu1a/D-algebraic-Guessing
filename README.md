# DalgGuessing

[![arXiv](https://img.shields.io/badge/arXiv-2510.26869-red.svg)](https://arxiv.org/abs/2510.26869)

> «The result of the mathematician's creative work is demonstrative reasoning, a proof, but the proof is discovered by plausible reasoning, by GUESSING.» — George Pólya

![Analogy](Analogy.png)
*The young fisherman caught a 'necklace' of consecutive terms of a sequence, providing strong evidence that this sequence belongs to a class of objects defined by its 'fishing rod' (being D-algebraic).*

## Introduction

This GitHub page is devoted to the Maple package **DalgGuessing**.

DalgGuessing is a standalone package for the guessing software from the Maple package [NLDE](https://github.com/T3gu1a/D-algebraic-functions) dedicated to operations with D-algebraic funtions and sequences.
The package currently exports two commands for recovering D-algebraic functions from the first few coefficients of their power series. These are:

- ``DalgFunGuess`` : guess an algebraic differential equation in characteristic zero from a list of initial coefficients of the D-algebraic function sought.
- ``modDalgFunGuess``: guess an algebraic differential equation over $\mathbb{Z} /m \mathbb{Z}$ from a list of initial coefficients of the D-algebraic function sought. Here $m$ is a positive integer (usually prime) chosen by the user.


## Requirements

A recent version of [Maple](https://www.maplesoft.com/) (version $\geq$ 2019).

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

- [Bertrand Teguia Tabuguia](https://bertrandteguia.com), [Email Me](mailto:email@bertrandteguia.com)
