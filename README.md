# SIF

A simple python module for conveniently building, modifying and converting LAMMPS input files.

This module was primarily made to ease the process of working with the REACTER module, but may have tools that can be useful for simulation setup in general.

## Overview

The main capabilities of SIF include:

- Simple iteration over atoms and topology types (bonds, angles, etc.)
- Tidy insertion and deletion of atoms, topologies and their types
- Assignment and management of type labels (introduced in LAMMPS 15Sep2022)
- Reading & writing LAMMPS .data files
- Reading & writing LAMMPS molecule files (used in e.g. the REACTER module)
- Tools for REACTER fragment file generation

## How-To

### Getting Started

To get up and running, add `import SIF` to the top of your python file.

A simulation system is referred to as a `World`. A World can be created from scratch using `SIF.World()`, or read from a file using functions available in the submodule `SIF.io`.

```python
import SIF
world_from_scratch = SIF.World(xlo=0.0, ylo=0.0, zlo=0.0, xhi=10.0, yhi=10.0, zhi=10.0)

import SIF.io
world_from_file = SIF.io.read_lammps_data("system.data")
```

From here, the World's atoms and topologies are available as well as their types:

- `World.atoms`, `World.atom_types`
- `World.bonds`, `World.bond_types`
- `World.angles`, `World.angle_types`
- `World.dihedrals`, `World.dihedral_types`
- `World.impropers`, `World.improper_types`(improper dihedrals)

### Adding & Removing Types

While you *can* edit the above arrays inplace if you want, this will tend to cause mismatches in type assignments etc, especially when changing the type arrays. Instead, functions are available to perform these operations while handling any side-effects correctly. For example, inserting an atom type will increase the atom type IDs of all those above the inserted type by 1. These functions are all accessible as `World.add_XXX` and `World.delete_XXX`. `atom`, `bond`, etc. and `atom_type`, `bond_type`, etc. are all valid in place of XXX.

### Automatic Type Labels

**NOTE:** This feature as-written is quite narrow. It works only on files produced by Moltemplate and currently applied only to the OPLS-AA forcefield (although it shouldn't be too difficult to apply to others).

### Examples

Example scripts showing off some of the potential uses of the module are in examples/ .  

- frag_gen/ : REACTER fragment/template file generation.
- topo_extension/ : Merging the topology lists of two separate input files.

## Installation

The module can be installed in either of two ways:

1) `pip install 'git+https://github.com/k-harris27/SIF'`

or

2)  i) `git clone https://github.com/k-harris27/SIF`  
   ii) `pip install ./SIF`
