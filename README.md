# BDMM-Flow

BDMM-Flow is a [BEAST 2](http://www.beast2.org/) package for
performing phylodynamic inference under multi-type birth-death models.

The package is heavily inspired by and designed to be a drop-in-replacement of the [BDMM-Prime](https://github.com/tgvaughan/BDMM-Prime) package. BDMM-Flow implements the same model, but uses a flow-based representation of the model likelihood. This representation is faster for large trees and was introduced in the following paper:

* Stilianos Louca, Matthew W Pennell, A General and Efficient Algorithm for the Likelihood of Diversification and Discrete-Trait Evolutionary Models, Systematic Biology, Volume 69, Issue 3, May 2020, Pages 545–556, [doi.org/10.1093/sysbio/syz055](https://doi.org/10.1093/sysbio/syz055)
§
This package adds additional modifications to the algorithm in order to further improve the performance.

## 🌴 Using the Package

In order to try out BDMM-Flow, add https://raw.githubusercontent.com/tochsner/BDMM-Flow/main/version.xml as a third party
BEAST package repository and installing the package that appears (see [this article](https://www.beast2.org/managing-packages/) for more information).

## 🔧 Building from Source

To build BDMM-Flow from source, you need to install the following dependencies:

- OpenJDK version 17 or greater
- the Apache Ant build system

Once these are installed and in your execution path, you can build the package by simply running the following command:

```sh
ant
```

## 👋 Acknowledgements

This package is heavily inspired by the [BDMM-Prime](https://github.com/tgvaughan/BDMM-Prime) package developed by [Tim Vaughan](https://github.com/tgvaughan). BDMM-Prime in turn is a fork of the original [BDMM package](https://github.com/denisekuehnert/bdmm) by [Denise Kühnert](https://github.com/denisekuehnert/)
and [Jérémie Scire](https://github.com/jscire).

The flow algorithm was developed by Louca and Pennell in this paper:

* Stilianos Louca, Matthew W Pennell, A General and Efficient Algorithm for the Likelihood of Diversification and Discrete-Trait Evolutionary Models, Systematic Biology, Volume 69, Issue 3, May 2020, Pages 545–556, [doi.org/10.1093/sysbio/syz055](https://doi.org/10.1093/sysbio/syz055)

## 📄 License

BDMM-Flow is free software.  It is distributed under the terms of version 3 of the GNU General Public License.  A copy of this license should be found in the file COPYING located in the root directory of this repository. If this file is absent for some reason, it can also be retrieved from https://www.gnu.org/licenses.