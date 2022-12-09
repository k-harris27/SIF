from utils import io,system,typing_tools

# Read in topology & settings (Forcefield) info of world to have atom/topo types added.
world = io.read_lammps_data("test_data/bulk_400D200M_s2.data","test_data/monomer.in.settings")

# Read in topology & settings (Forcefield) info of world used to take atom/topo types from (reference).
ref = io.read_lammps_data("test_data/poly_atoms.data","test_data/polymerisable.in.settings")

# Atom types have to be added manually currently...
for t in ((4,15.999),(5,1.008),(8,14.007),(9,14.007),(11,1.008)):
    world.add_atom_type(type_id=t[0],mass=t[1])

# Just checks to see the new atom types look right.
#   Prints "whoops" if there are any atoms of the new types
#   since we didn't make any actual atoms of those types.
print("Atoms")
print(world.atom_type_params)
[print("Whoops!") for a in world.atoms if a.type_id in (4,5,8,9,11)]  # Check atom types have been moved correctly
print("")

# Update topology types of world to match ref, then checks that there are no actual bonds of the new types.
# The print lines will break if the topology list breaks - things seem to work fine so can probably
#   delete them.
typing_tools.update_types_to_match(world,ref, debug = True)
[print(f"Whoops! (bond edition)") for b in world.bonds if b.type_id in (0,1)]
[print(f"Whoops! (angle edition)") for a in world.angles if a.type_id in (0,1,4,8,16)]
[print(f"Whoops! (dihedral edition)") for d in world.dihedrals if d.type_id in (0,1,2,3,4,5,9,15,16,17)]

# Finally, write type-extended data out
io.write_lammps_data(world, "test_data/bulk_400D200M_s2_ext.data")