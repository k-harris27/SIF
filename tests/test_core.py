import pytest
import SIF
from copy import deepcopy

@pytest.fixture(scope="function")
def empty_world() -> SIF.World:
    return SIF.World(-10, -10, -10, 10, 10, 10)

@pytest.fixture(scope="function")
def three_atom_world(empty_world : SIF.World) -> SIF.World:
    three_atoms = empty_world
    three_atoms.add_atom_type(12.0, "C")
    [three_atoms.add_atom(0, 0, (r, r, r)) for r in range(3)]
    return three_atoms

@pytest.fixture(scope="function")
def three_bonded_world(three_atom_world : SIF.World) -> SIF.World:
    world = three_atom_world
    for i in range(3):
        b_type = f"b{i}"
        world.add_bond_type(0.0, 0.0, type_name=f"b{i}")
        world.add_bond(i, (i+1) % 3, b_type)
    return world

def test_add_atom_type(empty_world : SIF.World):
    empty_world.add_atom_type(12.0, "C")
    assert empty_world.atom_types[0].name == "C"

def test_width(empty_world : SIF.World):
    assert all(w == 20 for w in empty_world.width)

def test_merge(three_atom_world : SIF.World):
    world1 = deepcopy(three_atom_world)
    world2 = three_atom_world
    world1.merge(world2, type_merge_style="Keep")
    assert world1.n_atoms == 6

def test_merge_with_type_merge(three_atom_world : SIF.World, empty_world : SIF.World):
    world1 = three_atom_world
    world1.add_atom_type(14.0, "N")  # World 1 has C & N
    world2 = empty_world
    world2.add_atom_type(12.0, "C")
    world2.add_atom_type(16.0, "O")  # World 2 has C & O
    world1.merge(world2, type_merge_style="Merge")
    assert [t.name for t in world1.atom_types] == ["C", "N", "O"]

def test_delete_bonded_atom(three_bonded_world : SIF.World):
    three_bonded_world.delete_atom(1)
    bonds = three_bonded_world.bonds
    bond_type_name = three_bonded_world.bond_types[bonds[0].type_id].name
    assert len(bonds) == 1 and  bond_type_name == "b2"

def test_get_bond_index(three_bonded_world : SIF.World):
    assert three_bonded_world.get_bond_index(1,2) == 1

def test_validate_topo_atoms(three_atom_world : SIF.World):
    assert three_atom_world._validate_topo_atoms(2,0,1) == [1,0,2]
