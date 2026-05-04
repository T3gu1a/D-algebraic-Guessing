# DalgGuessing

[![arXiv](https://img.shields.io/badge/arXiv-2510.26869-red.svg)](https://arxiv.org/abs/2510.26869)

> «The result of the mathematician's creative work is demonstrative reasoning, a proof, but the proof is discovered by plausible reasoning, by GUESSING.» — George Pólya

![Analogy](/images/Analogy.png)
*The young fisherman caught a 'necklace' of consecutive terms of a sequence, providing strong evidence that this sequence belongs to a class of objects defined by its 'fishing rod' (being D-algebraic).*

## Introduction

This GitHub page is devoted to the Maple package **DalgGuessing**.

DalgGuessing is a standalone package for the guessing software from the Maple package [NLDE](https://github.com/T3gu1a/D-algebraic-functions) dedicated to operations with D-algebraic funtions and sequences.
The package currently exports two commands for recovering D-algebraic functions from the first few coefficients of their power series. These are:

- ``DalgFunGuess`` : guess an algebraic differential equation in characteristic zero from a list of initial coefficients of the D-algebraic function sought.
- ``modDalgFunGuess``: guess an algebraic differential equation over $\mathbb{Z} /m \mathbb{Z}$ from a list of initial coefficients of the D-algebraic function sought. Here $m$ is a positive integer (usually prime) chosen by the user.


## Requirements

A recent version of [Maple](https://www.maplesoft.com/) (version $\geq$ 2019), best for versions $\geq 2023$.

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

## Calling sequence and arguments

``DalgFunGuess(L,degADE,degPoly,devars,startfromord,linsolver,inputConstants,allPolyDeg,termsOfDegPoly,sparsity,maxIteration)``

- ``L``: the list of initial terms from index $0$. This is the only mandatory arguments. The other arguments have default values.
- ``degADE``: the (total) degree of the ADE sought, with a default value of $2$.
- ``degPoly``: the maximum degree of the polynomial coefficients, with a default value of $2$.
- ``devars``: the dependent variable, with a default value of $y(x)$.
- ``startfromord``: the starting order of the search, with a default value of $0$.
- ``linsolver``: the method use to solve the underlying linear system with default value AlgebraicFunction for Maple's `SolveTools:-Linear`. The option `HardSystem` is an alternative which is more effective for larger systems or systems with many symbolic variables.
- ``inputConstants``: a set of parameters or symbolic variables appearing in the initial terms in `L`. Default value `{}`.
- ``allPolyDeg``: a boolean variable with a default value of `false`. When set to `true`, it prevents the code from stopping when it reaches the lowest order where the data in `L` becomes insufficient. Instead, it fixes that order and continues the search by allowing polynomial coefficients with degrees strictly less than `degPoly` to reduce the number of unknowns.
- ``termsOfDegPoly``: a positive integer representing the number of polynomial coefficients of degree `degPoly` to include if the data in `L` is insufficient for the general search performed when `allPolyDeg=true`. This can be used to reduce the order enough to allow the search to proceed. This is particularly useful when `allPolyDeg=false` fails to initiate any computation.
- ``sparsity``: a positive fraction between 0 (default) and 1. When non-zero with `allPolyDeg=true`, it tells the code to search for a differential equation of degree at most `degPoly` where a `sparsity` fraction of all polynomial coefficients are forced to be 0. This is the "desperate search" approach!
- ``maxIteration``: a positive integer or `infinity`. When finite, it bounds the number of iterations and randomizes the search when `allPolyDeg=true`.

Except for the `linsolver` argument, `modDalgFunGuess` shares the same arguments as `DalgFunGuess`, with the addition of:

- ``modulus``: a positive integer representing the modulus for the computation over $\mathbb{Z}/n\mathbb{Z}$, where $n$ is the `modulus`. Default value is $7$.

## Author

- [Bertrand Teguia Tabuguia](https://bertrandteguia.com), [Email Me](mailto:email@bertrandteguia.com)
