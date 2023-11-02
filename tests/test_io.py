import pytest
import SIF
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
    assert oplsaa_world.bond_types[0].name == "HO-OH"

def test_no_dupes_in_inferred(oplsaa_world : SIF.World):
    world = oplsaa_world
    avail_kinds = world._available_topo_types
    has_duplicates = {kind:None for kind in avail_kinds}
    for topo_kind in world._available_topo_types:
        topo_types = world._get_topo_type_list(topo_kind)
        type_names = [t.name for t in topo_types]
        # Fun trick - turn into a set (no duplicates).
        # If the lengths are different, we had duplicates.
        dupes = (len(type_names) != len(set(type_names)))
        has_duplicates[topo_kind] = dupes
        if dupes:
            print(f"Duplicate names: {[t for t in type_names if type_names.count(t) > 1]}")

    assert all([dupes is False for dupes in has_duplicates.values()])