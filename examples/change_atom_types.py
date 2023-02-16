from SIF import io, typing_tools
from SIF.system import lammps_index

world = io.read_lammps_data("test_data/dgeba-dm_frag.data")

typing_tools.change_atom_type(world,
                                new_type=lammps_index(13),
                                target_type=lammps_index(2),
                                target_indices=list(range(0,world.n_atoms,2)),
                                new_mass=12.011)

io.write_lammps_data(world,"atom_types_changed.data")