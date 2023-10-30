import SIF

world = SIF.io.read_lammps_data("test_data/dgeba-dm_frag.data")

world.add_atom_type(12.001, name="CE")

SIF.typing_tools.change_atom_type(world,
                                new_type=SIF.lammps_index(13),
                                target_type=SIF.lammps_index(2),
                                target_indices=list(range(0,world.n_atoms,2))
                                )

SIF.io.write_lammps_data(world,"test_data/OUT_atom_types_changed.data")