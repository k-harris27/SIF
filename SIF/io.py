from distutils.log import warn
from email import header
from .system import *
import resource

debug = False

def read_lammps_data(path_to_file : str, path_to_params : str = None) -> World:
    """Read LAMMPS .data file into a world object.\n
    path_to_file would be a .data file, path_to_params is optional and should direct to the file specifying interaction coefficients."""
    world = World()  # Create world instance to be populated

    num_atoms = None
    num_bonds = None
    num_angles = None
    num_dihedrals = None
    num_impropers = None

    warned_sections = []  # List of section names that unexpected warnings have been printed for. 

    section = "Headers"  # We start with headers
    with open(path_to_file,"r") as file:
        comment_line = True  # First line is comment
        for line in file:
            if comment_line:
                comment_line = False
                continue

            line,_,comment = line.partition("#")  # Separate code and comments
            line = line.strip()  # Remove all whitespace (AND NEWLINE CHARACTERS)

            if not line: continue  # Line is empty, skip.

            # Only section titles start with uppercase. Use title to indicate current section.
            if line[0].isupper():
                section = line
                continue

            
            words = line.split()

            if section == "Headers":
                
                if len(words) == 2: keyword = words[1]
                elif len(words) == 3: keyword = f"{words[1]} {words[2]}"
                elif len(words) == 4: keyword = f"{words[2]} {words[3]}"
                else: raise ValueError(f"Unexpected line in header: {line}")
    
                if keyword == "atoms":
                    num_atoms = int(words[0])
                    world.atoms = [None] * num_atoms
                elif keyword == "bonds":
                    num_bonds = int(words[0])
                    world.bonds = [None] * num_bonds
                elif keyword == "angles":
                    num_angles = int(words[0])
                    world.angles = [None] * num_angles
                elif keyword == "dihedrals":
                    num_dihedrals = int(words[0])
                    world.dihedrals = [None] * num_dihedrals
                elif keyword == "impropers":
                    num_impropers = int(words[0])
                    world.impropers = [None] * num_impropers
                elif keyword == "atom types":
                    world.append_atom_types([0.] * int(words[0]))
                elif keyword == "bond types":
                    world.append_bond_types([None] * int(words[0]))
                elif keyword == "angle types":
                    world.append_angle_types([None] * int(words[0]))
                elif keyword == "dihedral types":
                    world.append_dihedral_types([None] * int(words[0]))
                elif keyword == "improper types":
                    world.append_improper_types([None] * int(words[0]))
                elif keyword == "xlo xhi":
                    world.set_dims(xlo = float(words[0]), xhi = float(words[1]))
                elif keyword == "ylo yhi":
                    world.set_dims(ylo = float(words[0]), yhi = float(words[1]))
                elif keyword == "zlo zhi":
                    world.set_dims(zlo = float(words[0]), zhi = float(words[1]))
                else:
                    raise ValueError(f"Unexpected line in header: {line}")

            elif section == "Masses":
                world.set_atom_type_mass(type_index = int(words[0])-1,  # LAMMPS indices start at 1
                mass = float(words[1]))

            elif section == "Atoms":
                """
                Atoms section words:
                0 = atom ID
                1 = mol ID (generally useless?)
                2 = atom type
                3 = charge
                4,5,6 = x,y,z position
                7,8,9 = image flags
                """
                atom_id = int(words[0])-1
                atom_type = int(words[2])-1
                charge = float(words[3])
                pos = [float(r) for r in words[4:7]]
                img_flags = [int(f) for f in words[7:10]]
                if len(img_flags) == 3:
                    pos = [p + img*world.width[dim] for p,img,dim in zip(pos,img_flags,range(3))]
                world.atoms[atom_id] = Atom(atom_type,charge,pos)

            elif section == "Bonds":
                """
                Bonds section words:
                0 = Bond number ID
                1 = Bond type
                2 = atom 1 ID
                3 = atom 2 ID
                """
                i = int(words[0])-1
                bond_type = int(words[1])-1
                atom_1 = int(words[2])-1
                atom_2 = int(words[3])-1
                atom_1,atom_2 = world._validate_topo_atoms(atom_1,atom_2)
                world.bonds[i] = Topology(atom_1, atom_2, type_id = bond_type)

            elif section == "Angles":
                """
                Angles section words:
                0 = Angle number ID
                1 = Angle type
                2,3,4 = atom IDs
                """
                i = int(words[0])-1
                angle_type = int(words[1])-1
                atoms = world._validate_topo_atoms(*(int(w)-1 for w in words[2:]))
                world.angles[i] = Topology(*atoms, type_id = angle_type)

            elif section == "Dihedrals":
                """
                Dihedrals section words:
                0 = Dih number ID
                1 = Dih type
                2,3,4,5 = atom IDs
                """
                i = int(words[0])-1
                dih_type = int(words[1])-1
                atoms = world._validate_topo_atoms(*(int(w)-1 for w in words[2:]))
                world.dihedrals[i] = Topology(*atoms, type_id = dih_type)

            elif section == "Impropers":
                """
                Impropers section words:
                0 = Imp number ID
                1 = Imp type
                2,3,4,5 = atom IDs
                """
                i = int(words[0])-1
                imp_type = int(words[1])-1
                atoms = world._validate_topo_atoms(*(int(w)-1 for w in words[2:]))
                world.impropers[i] = Topology(*atoms, type_id = imp_type)

            elif section not in warned_sections:
                print(f"WARNING: Unexpected section title in LAMMPS input file: {section}")
                warned_sections.append(section)

    # Interpret coefficients file, if it exists
    if path_to_params:
        with open(path_to_params,"r") as file:
            for line in file:
                line,_,comment = line.partition("#")  # Separate code and comments
                line = line.strip()  # Remove all whitespace (AND NEWLINE CHARACTERS)

                if not line: continue  # Line is empty, skip.

                words = line.split()

                if words[0] == "pair_coeff":
                    pass  # TODO: Storing pair coefficients???
                # Collect and store topology parameters. There's better ways to do this but eh.
                for topo_name in ("bond","angle","dihedral","improper"):
                    if words[0] == f"{topo_name}_coeff":
                        topo_params = world._get_topo_param_list(topo_name)
                        id = int(words[1])-1
                        params = [float(w) for w in words[2:]]
                        topo_params[id] = params
                

    if debug: _debug_print_resource()
    return world

def write_lammps_data(world : World, path_to_file : str, comment : str = "") -> None:
    """Write world to file in LAMMPS .data format"""

    with open(path_to_file, "w") as file:
        # Comments
        file.write(comment + "\n")
        file.write("\n")

        # Header section
        file.write(f"{len(world.atoms)} atoms\n")
        file.write(f"{len(world.bonds)} bonds\n")
        file.write(f"{len(world.angles)} angles\n")
        file.write(f"{len(world.dihedrals)} dihedrals\n")
        file.write(f"{len(world.impropers)} impropers\n")
        file.write("\n")
        file.write(f"{len(world.atom_type_params)} atom types\n")
        file.write(f"{len(world.bond_params)} bond types\n")
        file.write(f"{len(world.angle_params)} angle types\n")
        file.write(f"{len(world.dihedral_params)} dihedral types\n")
        file.write(f"{len(world.improper_params)} improper types\n")
        file.write("\n")
        file.write(f"{world.xlo} {world.xhi} xlo xhi\n")
        file.write(f"{world.ylo} {world.yhi} ylo yhi\n")
        file.write(f"{world.zlo} {world.zhi} zlo zhi\n")
        file.write("\n")

        # Masses section
        file.write("Masses\n")
        file.write("\n")
        for i,m in enumerate(world.atom_type_params):
            file.write(f"{i+1} {m}\n")
        file.write("\n")

        # Atoms section
        file.write("Atoms\n")
        file.write("\n")
        for i,a in enumerate(world.atoms):
            file.write(f"{i+1} {1} {a.type_id+1} {a.charge} {' '.join([str(round(r,6)) for r in a.pos])}\n")
        file.write("\n")

        # Topology sections
        for topo_name in ("bond","angle","dihedral","improper"):
            topo = world._get_topo_list(topo_name)
            header_text = topo_name.capitalize() + "s"
            if len(topo) > 0: _write_lammps_topo(header_text,topo,file)


def read_react_template(path_to_file : str) -> World:
    world = World()  # Create world instance to be populated

    section = "Headers"  # We start with headers
    with open(path_to_file,"r") as file:
        comment_line = True  # First line is comment
        n_types = {"atom":0,"bond":0,"angle":0,"dihedral":0,"improper":0}
        for line in file:
            if comment_line:
                comment_line = False
                continue

            line,_,comment = line.partition("#")  # Separate code and comments
            line = line.strip()  # Remove all whitespace (AND NEWLINE CHARACTERS)

            if not line: continue  # Line is empty, skip.

            # Only section titles start with uppercase. Use title to indicate current section.
            if line[0].isupper():
                section = line
                continue

            words = line.split()

            if section == "Headers":
                val = int(words[0])
                if words[1] == "atoms":
                    world.atoms = [Atom(0,0,[0,0,0]) for _ in range(val)]
                elif words[1] in ("bonds", "angles", "dihedrals", "impropers"):
                    topo_name = words[1][:-1]
                    topo = world._get_topo_list(topo_name)
                    topo.extend([None] * val)  # Uses mutability to set list (starts empty)
                else:
                    raise ValueError(f"Unexpected line in header: {line}")
            elif section == "Coords":
                i = int(words[0])-1
                pos = [float(w) for w in words[1:]]
                world.atoms[i].pos = pos
            elif section == "Types":
                i = int(words[0])-1
                type_id = int(words[1])-1
                n_types["atom"] = max(n_types["atom"],type_id+1)
                world.atoms[i].type_id = int(words[1])-1
            elif section == "Charges":
                i = int(words[0])-1
                world.atoms[i].charge = float(words[1])
            elif section == "Fragments":
                frag_id = words[0]
                atom_ids = [int(w)-1 for w in words[1:]]
                world.add_fragment(frag_id,atom_ids)
            elif section in ("Bonds", "Angles", "Dihedrals", "Impropers"):
                topo_name = section.lower()[:-1]
                topo = world._get_topo_list(topo_name)
                i = int(words[0])-1
                topo_type = int(words[1])-1
                n_types[topo_name] = max(n_types[topo_name],topo_type+1)
                atoms = world._validate_topo_atoms(*(int(w)-1 for w in words[2:]))
                topo[i] = Topology(*atoms, type_id = topo_type)
            else:
                raise ValueError(f"Unexpected section title in react template file: {section}")

            world.atom_type_params = [1.0] * n_types["atom"]

        for topo_name in ("bond","angle","dihedral","improper"):
            param_list = world._get_topo_param_list(topo_name)
            print(f"{topo_name}: {n_types[topo_name]}")
            param_list.extend([None] * n_types[topo_name])
    
    return world

def write_react_template(world : World, path_to_file : str, comment : str = "") -> None:
    """
    Write world to file in react .template format
    """

    with open(path_to_file, "w") as file:
        # Comments
        file.write(comment + "\n")
        file.write("\n")

        # Header section
        file.write(f"{len(world.atoms)} atoms\n")
        file.write(f"{len(world.bonds)} bonds\n")
        file.write(f"{len(world.angles)} angles\n")
        file.write(f"{len(world.dihedrals)} dihedrals\n")
        file.write(f"{len(world.impropers)} impropers\n")
        if len(world.fragments) > 0: file.write(f"{len(world.fragments)} fragments\n")
        file.write("\n")

        # Coords section
        file.write("Coords\n")
        file.write("\n")
        for i,a in enumerate(world.atoms):
            file.write(f"{i+1} {' '.join([str(round(r,6)) for r in a.pos])}\n")
        file.write("\n")

        # Types section
        file.write("Types\n")
        file.write("\n")
        for i,a in enumerate(world.atoms):
            file.write(f"{i+1} {a.type_id+1}\n")
        file.write("\n")

        # Charges section
        file.write("Charges\n")
        file.write("\n")
        for i, a in enumerate(world.atoms):
            file.write(f"{i+1} {a.charge}\n")
        file.write("\n")

        # Fragments section
        if len(world.fragments) > 0:
            file.write("Fragments\n")
            file.write("\n")
            for frag_id,atom_ids in world.fragments.items():
                lammps_ids = [str(i+1) for i in atom_ids]
                file.write(f"{frag_id} {' '.join(lammps_ids)}\n")
            file.write("\n")

        # Topology sections
        for topo_name in ("bond","angle","dihedral","improper"):
            topo = world._get_topo_list(topo_name)
            header_text = topo_name.capitalize() + "s"
            if len(topo) > 0: _write_lammps_topo(header_text,topo,file)

        
# Used to write any topo sections for lammps-style files
def _write_lammps_topo(topo_name, topo, file):
    file.write(topo_name + "\n")
    file.write("\n")
    for i,t in enumerate(topo):
        file.write(f"{i+1} {t.type_id+1} {' '.join([str(a+1) for a in t])}\n")
    file.write("\n")

def _debug_print_resource():
    print('Peak Memory Usage =', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    print('User Mode Time =', resource.getrusage(resource.RUSAGE_SELF).ru_utime)
    print('System Mode Time =', resource.getrusage(resource.RUSAGE_SELF).ru_stime)
