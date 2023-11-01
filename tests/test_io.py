import pytest
import SIF
import SIF.io
from SIF.forcefields import ForceFields
from pathlib import Path

def test_forcefield():
    dir = f"{Path(__file__).resolve().parents[1]}/test_data/poly_atoms.data"
    world = SIF.io.read_lammps_data(dir, type_name_lookup=ForceFields.oplsaa())
    available_types = ForceFields.oplsaa().values()
    world_types = [t.name for t in world.atom_types]
    assert all(t in available_types for t in world_types)