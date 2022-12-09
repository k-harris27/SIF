from utils import system
from utils import io
from utils import template_tools
from copy import copy

"""
Reading lammps data
world = io.read_lammps_data("test_data/system.data")

[print(a.charge) for a in world.atoms[:5]]

print(world.bonds[:5])

print(world.angles[1200:1205])

print(world.dihedrals[:5])

print(world.impropers[:5])
"""

"""
#Merging Template Files

dgeba = io.read_lammps_data("test_data/dgeba_frag.data")

mxda = io.read_react_template("test_data/mxda_frag.template")

dgeba.merge(mxda)

io.write_react_template(dgeba, "test_data/dgeba-mxda_frag.template", "DGEBA MXDA Merged Template")
io.write_lammps_data(dgeba, "test_data/dgeba-mxda_frag.data", "DGEBA MXDA Merged Template")

dgeba = io.read_lammps_data("test_data/dgeba_frag.data")

dm = io.read_lammps_data("test_data/dm_frag.data")

dgeba.merge(dm)

io.write_react_template(dgeba, "test_data/dgeba-dm_frag.template", "DGEBA DM-dimer Merged Template")
io.write_lammps_data(dgeba, "test_data/dgeba-dm_frag.data", "DGEBA DM-dimer Merged Template")
"""


#Test template creation from full data

world = io.read_lammps_data("test_data/system.data")

template_tools.get_connected_atoms(world, 0)


"""
data to template

dgeba = io.read_lammps_data("test_data/dgeba_frag.data")
io.write_react_template(dgeba, "test_data/dgeba_frag.template", comment="DGEBA Template with fixed positions")
dm = io.read_lammps_data("test_data/dm_frag.data")
io.write_react_template(dm, "test_data/dm_frag.template", comment="DM Template with fixed positions")
"""

