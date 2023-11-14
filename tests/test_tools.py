import pytest
import SIF
from SIF.template_tools import get_connected_atoms
from SIF.typing_tools import update_types_to_match
from SIF.forcefields import ForceFields
from pathlib import Path
from copy import deepcopy

@pytest.fixture
def oplsaa_world():
    dir = f"{Path(__file__).resolve().parent}/test_data/oplsaa_atoms.data"
    return SIF.io.read_lammps_data(dir, forcefield=ForceFields.oplsaa)

def test_bonded_search_1_atom(oplsaa_world : SIF.World):
    world = oplsaa_world
    sel = get_connected_atoms(world, 21, extent=3)
    assert sel.n_atoms == 10  # 1 DGEBA reactive fragment

def test_bonded_search_2_atoms(oplsaa_world : SIF.World):
    world = oplsaa_world
    sel = get_connected_atoms(world, (21,49), extent=3)
    assert sel.n_atoms == 19  # DGEBA + MXDA reactive fragments

def test_update_atom_types_to_match(oplsaa_world : SIF.World):
    target = oplsaa_world
    ref = deepcopy(oplsaa_world)

    # Delete all atoms and topology so we don't get anything
    #   oointing to types that don't exist.
    target.atoms = []
    for topo_kind in target._available_topo_types:
        target._get_topo_list(topo_kind).clear()

    for i in (1,6,7):
        target.atom_types.pop(i)

    update_types_to_match(target, ref)

    target_names = [t.name for t in target.atom_types]
    ref_names = [t.name for t in ref.atom_types]
    assert target_names == ref_names

def test_update_topo_types_to_match(oplsaa_world : SIF.World):
    target = oplsaa_world
    ref = deepcopy(oplsaa_world)

    # Delete all atoms and topology so we don't get anything
    #   oointing to types that don't exist.
    target.atoms = []
    for topo_kind in target._available_topo_types:
        target._get_topo_list(topo_kind).clear()

    for i in (1,6,7):
        target.dihedral_types.pop(i)
    
    update_types_to_match(target, ref)

    target_names = [t.name for t in target.dihedral_types]
    ref_names = [t.name for t in ref.dihedral_types]

    assert target_names == ref_names