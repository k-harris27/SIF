from ._core import *
import resource
from typing import Dict, Tuple
from SIF.forcefields import ForceField

# TODO: Read/Write molecule files properly too.

def read_lammps_data(path_to_file : str, path_to_params : str = None, forcefield : ForceField = None) -> World:
    """Read LAMMPS .data file into a world object.\n
    path_to_file would be a .data file, path_to_params is optional and should direct to the file specifying interaction coefficients."""
    world = World()  # Create world instance to be populated

    num_atoms = None
    num_bonds = None
    num_angles = None
    num_dihedrals = None
    num_impropers = None

    if forcefield is not None: type_name_lookup = forcefield.atom_names

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
                    world.append_atom_types([AtomType(mass=0) for _ in range(int(words[0]))])
                elif keyword == "bond types":
                    world.append_bond_types([TopologyType() for _ in range(int(words[0]))])
                elif keyword == "angle types":
                    world.append_angle_types([TopologyType() for _ in range(int(words[0]))])
                elif keyword == "dihedral types":
                    world.append_dihedral_types([TopologyType() for _ in range(int(words[0]))])
                elif keyword == "improper types":
                    world.append_improper_types([TopologyType() for _ in range(int(words[0]))])
                elif keyword == "xlo xhi":
                    world.set_dims(xlo = float(words[0]), xhi = float(words[1]))
                elif keyword == "ylo yhi":
                    world.set_dims(ylo = float(words[0]), yhi = float(words[1]))
                elif keyword == "zlo zhi":
                    world.set_dims(zlo = float(words[0]), zhi = float(words[1]))
                else:
                    raise ValueError(f"Unexpected line in header: {line}")

            elif "Type Labels" in section:
                if len(words) != 2:
                    raise ValueError(f"Invalid atom type label input: {line}")
                type_id = int(words[0])-1
                if section == "Atom Type Labels":
                    world.atom_types[type_id].name = words[1]
                    continue
                # Probably not very efficient but there shouldn't be many types anyway...
                for topo_kind in world._available_topo_types:
                    if section == f"{topo_kind.capitalize()} Type Labels":
                        topo_types = world._get_topo_type_list(topo_kind)
                        type_id = int(words[0])-1
                        topo_types[type_id].name = words[1]

            elif section == "Masses":
                world.atom_types[int(words[0])-1].mass = float(words[1])
                comment_words = comment.split()
                if len(comment_words) == 1:  # Assume it's an atom type string
                    name = comment_words[0]
                    if type_name_lookup is not None:
                        name = type_name_lookup[name]
                    world.atom_types[int(words[0])-1].name = name

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
                logger.warning(f"Unexpected section title in LAMMPS input file: {section}")
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
                        topo_types = world._get_topo_type_list(topo_name)
                        id = int(words[1])-1
                        params = [float(w) for w in words[2:]]
                        topo_types[id].parameters = params
                
    # Infer topology names if all atom type names are given.
    if all(t.name is not None for t in world.atom_types):
        world.infer_topo_names_from_atoms(override=False, forcefield=forcefield)

    _debug_print_resource()
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
        file.write(f"{len(world.atom_types)} atom types\n")
        file.write(f"{len(world.bond_types)} bond types\n")
        file.write(f"{len(world.angle_types)} angle types\n")
        file.write(f"{len(world.dihedral_types)} dihedral types\n")
        file.write(f"{len(world.improper_types)} improper types\n")
        file.write("\n")
        file.write(f"{world.xlo} {world.xhi} xlo xhi\n")
        file.write(f"{world.ylo} {world.yhi} ylo yhi\n")
        file.write(f"{world.zlo} {world.zhi} zlo zhi\n")
        file.write("\n")

        all_atom_names_defined = all(t.name is not None for t in world.atom_types)
        if all_atom_names_defined:
            file.write("Atom Type Labels\n\n")
            for i,t in enumerate(world.atom_types):
                file.write(f"{i+1} {t.name}\n")
            file.write("\n")
        
            # Override = False will prevent this from overwriting
            # if we've assigned manually/previously.
            world.infer_topo_names_from_atoms(override=False)
            for topo_kind in world._available_topo_types:
                file.write(f"{topo_kind.capitalize()} Type Labels\n\n")
                for i, topo_type in enumerate(world._get_topo_type_list(topo_kind)):
                    file.write(f"{i+1} {topo_type.name}\n")
                file.write("\n")

        # Masses section
        file.write("Masses\n")
        file.write("\n")
        for i,t in enumerate(world.atom_types):
            file.write(f"{i+1} {t.mass}\n")
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
        type_labels = {key:[] for key in n_types.keys()}
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
                    topo_kind = words[1][:-1]
                    topo = world._get_topo_list(topo_kind)
                    topo.extend([None] * val)  # Uses mutability to set list (starts empty)
                else:
                    raise ValueError(f"Unexpected line in header: {line}")
            elif section == "Coords":
                i = int(words[0])-1
                pos = [float(w) for w in words[1:]]
                world.atoms[i].pos = pos
            elif section == "Types":
                i = int(words[0])-1
                type_id, n_types["atom"] = _read_react_type(words[1], n_types["atom"], type_labels["atom"])
                world.atoms[i].type_id = type_id
            elif section == "Charges":
                i = int(words[0])-1
                world.atoms[i].charge = float(words[1])
            elif section == "Fragments":
                frag_id = words[0]
                atom_ids = [int(w)-1 for w in words[1:]]
                world.add_fragment(frag_id,atom_ids)
            elif section in ("Bonds", "Angles", "Dihedrals", "Impropers"):
                topo_kind = section.lower()[:-1]
                topo = world._get_topo_list(topo_kind)
                i = int(words[0])-1
                topo_type, n_types[topo_kind] = _read_react_type(words[1], n_types[topo_kind], type_labels[topo_kind])
                atoms = world._validate_topo_atoms(*(int(w)-1 for w in words[2:]))
                topo[i] = Topology(*atoms, type_id = topo_type)
            else:
                raise ValueError(f"Unexpected section title in react template file: {section}")

        world.atom_types = [AtomType(0.0, name) for name in type_labels["atom"]]

        for topo_kind in ("bond","angle","dihedral","improper"):
            param_list = world._get_topo_type_list(topo_kind)
            logger.debug(f"{topo_kind}: {n_types[topo_kind]}")
            param_list.extend(TopologyType(0.0, name=name) for name in type_labels[topo_kind])
    
    return world

def write_react_template(world : World, path_to_file : str, comment : str = "") -> None:
    """
    Write world to file in react .template format
    """

    all_type_labels_defined = all(t.name is not None for t in world.atom_types)
    if all_type_labels_defined:
        world.infer_topo_names_from_atoms(override=False)  # Override=False only updates unset type labels

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
            t = a.type_id + 1  # Conversion to LAMMPS type.
            if all_type_labels_defined:
                t =  world.atom_types[a.type_id].name  # Use atom type name if possible.
            file.write(f"{i+1} {t}\n")
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
        for topo_kind in world._available_topo_types:
            topo = world._get_topo_list(topo_kind)
            type_names = [t.name for t in world._get_topo_type_list(topo_kind)]
            header_text = topo_kind.capitalize() + "s"
            if len(topo) > 0: _write_lammps_topo(header_text,topo,file, type_names = type_names)


# Used to write any topo sections for lammps-style files
def _write_lammps_topo(header:str, topo:List[Topology], file, type_names = None):
    file.write(header + "\n")
    file.write("\n")
    for i,t in enumerate(topo):
        tp = t.type_id+1
        if type_names is not None:
            tp = type_names[t.type_id]
        file.write(f"{i+1} {tp} {' '.join([str(a+1) for a in t])}\n")
    file.write("\n")

def _read_react_type(read_type : str, n_types : int, type_labels : List[str]) -> Tuple[int,int]:
    """Works for atoms and topology. Takes in the type read from template file,
        and returns the new total number of types as well as the atom type number.
        
        Parameters
        ----------
        read_type : str
            Type read from file (may be a digit string if type labels aren't used).
        n_types : int
            Current minimum number of types.
        type_labels : list[str]
            Current list of all identified type labels.

        Returns
        -------
        type_id : int
            Input atom type index.
        new_n_types : int
            New minimum number of types.
    """
    
    if read_type.isdigit():
        type_id = int(read_type)-1
        n_types = max(n_types,type_id+1)
    elif read_type not in type_labels:
        type_id = n_types
        n_types += 1
        type_labels.append(read_type)
    else:
        type_id = type_labels.index(read_type)
    return type_id, n_types


def _debug_print_resource():
    logger.debug('Peak Memory Usage =', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    logger.debug('User Mode Time =', resource.getrusage(resource.RUSAGE_SELF).ru_utime)
    logger.debug('System Mode Time =', resource.getrusage(resource.RUSAGE_SELF).ru_stime)
