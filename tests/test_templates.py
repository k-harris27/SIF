import pytest
import SIF
from SIF.template_tools import get_connected_atoms
from SIF.forcefields import ForceFields
from pathlib import Path

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