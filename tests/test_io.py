import pytest
import SIF
import SIF.io
from SIF.forcefields import ForceFields
import os
from pathlib import Path

@pytest.fixture
def lmp_labelled_world():
    dir = f"{Path(__file__).resolve().parent}/test_data/type_labels.data"
    return SIF.io.read_lammps_data(dir)

@pytest.fixture
def out_world_path():
    fpath = f"{Path(__file__).resolve().parent}/test_data/__out.data"
    yield fpath  # Send the path to the test function, then...
    # Once test is done, delete test output file.
    os.remove(fpath)

@pytest.fixture
def oplsaa_world():
    dir = f"{Path(__file__).resolve().parent}/test_data/oplsaa_atoms.data"
    return SIF.io.read_lammps_data(dir, type_name_lookup=ForceFields.oplsaa())

def test_read_labelled_lammps_data(lmp_labelled_world : SIF.World):
    world = lmp_labelled_world
    assert world.atom_types[0].name == "A" and world.atom_types[1].name == "B" and world.bond_types[0].name == "A-B"

def test_write_labelled_lammps_data(lmp_labelled_world : SIF.World, out_world_path : str):
    in_world = lmp_labelled_world
    in_world_atom_names = [t.name for t in in_world.atom_types]
    SIF.io.write_lammps_data(in_world, out_world_path)
    reread_world = SIF.io.read_lammps_data(out_world_path)
    reread_world_atom_names = [t.name for t in in_world.atom_types]
    assert reread_world_atom_names == in_world_atom_names

def test_ff_atom_names(oplsaa_world : SIF.World):
    available_types = ForceFields.oplsaa().values()
    world_types = [t.name for t in oplsaa_world.atom_types]
    assert all(t in available_types for t in world_types)

def test_ff_bond_names(oplsaa_world : SIF.World):
    assert oplsaa_world.bond_types[0].name == "CA-HA"