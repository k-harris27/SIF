### Defining classes for system description ###

# TODO: Run test script creating/deleting etc??

# NOTE: Do I need this???
from argparse import ArgumentError
from multiprocessing.sharedctypes import Value
from typing import Iterable,List,Dict
import math
from copy import deepcopy


class Vector3:

    # Initialiser
    # If creating a vector3 from a list, use *list as single arg
    def __init__(self,x,y,z):
        self._vec = [x,y,z]

    def __repr__(self):
        return '{self.__class__.__name__}(x={self._vec[0]}, y={self._vec[1]}, z={self._vec[2]})'.format(self=self)


# TODO: Probably switch to getter & setter properties?
# https://stackoverflow.com/questions/2627002/whats-the-pythonic-way-to-use-getters-and-setters

class Atom:
    """
    A Class representing an individual atom

    Attributes
    ----------
    type_id : int
        The atom type index
    charge : float
        Atom partial charge
    pos : list of float
        List of length 3, defining x, y and z position of atom respectively
    img_flags : list of int
        List of length 3, defining x, y and z image flags of atom respectively
    """

    # Initialiser
    # type_id = lammps atom type
    # charge = atom partial charge (float)
    # pos = x,y,z position (list/array of 3 floats)
    # img_flags = x,y,z image flags (list/array of 3 integers)
    def __init__(self, type_id, charge, pos, img_flags=[0,0,0]):

        # Check type_id is an integer
        try:
            type_id = int(type_id)
        except ValueError as err:
            raise Exception(f"Atom type ID given as {type_id}, should be an integer.") from err

        # Check charge is a float
        try:
            charge = float(charge)
        except ValueError as err:
            raise Exception(f"Atom charge given as {charge}, should be a float.") from err

        # Check pos is a list-type variable
        try:
            pos = list(pos)
        except ValueError as err:
            raise Exception(f"Atom position {pos} should be a list or tuple.") from err
        
        # Check pos has 3 elements
        try:
            assert len(pos) == 3
        except AssertionError as err:
            raise Exception(f"Atom position {pos} should have exactly 3 elements.") from err

        # Check pos is all floats
        try:
            pos = [float(f) for f in pos]
        except ValueError as err:
            raise Exception(f"Atom position {pos} contains a non-float element.") from err

        # Check img_flags is a list-type variable
        try:
            img_flags = list(img_flags)
        except ValueError as err:
            raise Exception(f"Atom image flags {img_flags} should be a list or tuple.") from err

        # Check img_flags has 3 elements
        try:
            assert len(img_flags) == 3
        except AssertionError as err:
            raise Exception(f"Atom image flags {img_flags} should have exactly 3 elements.") from err

        # Check img_flags is all ints
        try:
            img_flags = [int(i) for i in img_flags]
        except ValueError as err:
            raise Exception(f"Atom image flags {img_flags} contains a non-integer element.") from err


        # Attributes
        self.type_id = type_id
        self.charge = charge
        self.pos = pos
        self.img_flags = img_flags
    
    def __copy__(self):
        """Ensures pos and image flag lists are copied properly too."""
        return Atom(self.type_id,self.charge,self.pos[:],self.img_flags[:])


class Topology(list):
    """
    A class used to represent a topology instance (a bond, an angle, ...), inheriting from list.

    By inheriting list, we get to say Topology[i] to get atom i of the topology.

    Attributes
    ----------
    type_id : int
        The type index of the given topology (used to define interaction parameters)
    
    Methods
    -------
    get_atoms():
        Returns an ordered list of atoms defining the topology
    """
    def __init__(self, *atoms, type_id=0):
        self.type_id = int(type_id)
        super().__init__(atoms)
    
    def __str__(self):
        return f"Topo(Type: {self.type_id}, Atoms: {self.get_atoms()})"

    def __repr__(self):
        return f"Topology({', '.join(str(i) for i in self.get_atoms())}, type_id={self.type_id})"
        
    def get_atoms(self):
        """
        Return an ordered list of atoms defining the topology

        Returns
        -------
        List of Atom : Ordered list of atoms defining the topology
        """
        return list(self)


class World:
    
    def __init__(self, xlo=0., ylo=0., zlo=0., xhi=0., yhi=0., zhi=0.):
        
        # Box Dimensions
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi

        # Atoms & topology
        self.atoms = []
        self.bonds = []  
        self.angles = []
        self.dihedrals = []
        self.impropers = []

        # Force field data
        self.atom_type_params = []  # atom_types contains mass of atom type i
        self.bond_params = []
        self.angle_params = []
        self.dihedral_params = []
        self.improper_params = []

        # Additional metadata
        self.fragments = {}  # Dict[str,list[int]]. Used by LAMMPS molecule files.

    @property
    def n_atom_types(self):
        """Return the number of atom types available in the system."""
        return len(self.atom_type_params)

    @property
    def width(self):
        return (self.xhi - self.xlo, self.yhi - self.ylo, self.zhi - self.zlo)

    @property
    def xwidth(self):
        return self.width[0]
    
    @property
    def ywidth(self):
        return self.width[1]
    
    @property
    def zwidth(self):
        return self.width[2]

    @property
    def total_charge(self):
        """Returns sum of charges of atoms in the World to 12 decimal places
        (rounded to remove floating point errors)"""
        return round(math.fsum([a.charge for a in self.atoms]),12)

    def add_fragment(self, frag_id : str = "", atom_ids : List[int] = []):
        """Assign a fragment grouping to the atoms with ids in atom_ids, giving the fragment the name frag_id.\n
        Fragments can be accessed from World.fragments.
        
        Arguments
        ---------
        frag_id : str
            The name to assign the fragment. If not given, assigns frag ID as first integer not in use.
        atom_ids : List[int]
            The list of atom IDs to include in the fragment.
        """

        # Check we have a string for the frag id.
        try:
            frag_id = str(frag_id)
        except:
            raise TypeError(f"frag_id must be a string or castable to a string. Received type {type(frag_id)}.")

        # If frag_id isn't given, find the first integer not in use.
        # Remove all leading/trailing whitespace, then check if string contains anything
        # (I.e. check if string doesn't contain any non-whitespace characters.)
        if not frag_id.strip():
            n = 1
            while True:
                frag_id = str(n)
                if frag_id in self.fragments.keys(): n+=1
                else: break  # current frag id (=str(n)) isn't in use, so keep it and exit the while loop.

        # Check we have a list of integers for atom_ids.
        assert isinstance(atom_ids,list), f"atom_ids must be a list or other iterable-type. Received type {type(atom_ids)}."
        for i in atom_ids:
            assert isinstance(i,int), f"elements in atom_ids must be integers. Received type {type(i)}."

        if frag_id in self.fragments.keys(): raise ArgumentError(f"A fragment with ID {frag_id} has already been assigned.")

        self.fragments[frag_id] = deepcopy(atom_ids)
    
    def set_dims(self, xlo=None, ylo=None, zlo=None, xhi=None, yhi=None, zhi=None):
        """Don't need to provide all, just set box dims of interest."""
        if xlo: self.xlo = xlo
        if xhi: self.xhi = xhi
        if ylo: self.ylo = ylo
        if yhi: self.yhi = yhi
        if zlo: self.zlo = zlo
        if zhi: self.zhi = zhi

    def merge(self, world_2) -> None:
        """
        Merge a second world into self, arranging topology and atoms automatically.
        """
        # Sort atom/topo types first
        atom_offset = len(self.atoms)

        n_params = [{   "atoms": self.n_atom_types,
                        "bonds": len(self.bond_params),
                        "angles": len(self.angle_params),
                        "dihedrals": len(self.dihedral_params),
                        "impropers": len(self.improper_params)},
                    {   "atoms": world_2.n_atom_types,
                        "bonds": len(world_2.bond_params),
                        "angles": len(world_2.angle_params),
                        "dihedrals": len(world_2.dihedral_params),
                        "impropers": len(world_2.improper_params)}]

        unequal_types = [n_params[0][t] != n_params[1][t] for t in ("atoms", "bonds", "angles", "dihedrals", "impropers")]
        if True in unequal_types: print("WARNING: Merged systems have different numbers of atom and/or topology types.")

        if n_params[1]["atoms"] > n_params[0]["atoms"]: self.atom_type_params = world_2.atom_type_masses.copy()
        if n_params[1]["bonds"] > n_params[0]["bonds"]: self.bond_params = world_2.bond_params.copy()
        if n_params[1]["angles"] > n_params[0]["angles"]: self.angle_params = world_2.angle_params.copy()
        if n_params[1]["dihedrals"] > n_params[0]["dihedrals"]: self.dihedral_params = world_2.dihedral_params.copy()
        if n_params[1]["impropers"] > n_params[0]["impropers"]: self.improper_params = world_2.improper_params.copy()

        self.atoms.extend(
            [Atom(a.type_id,a.charge,a.pos,a.img_flags.copy()) for a in world_2.atoms]
            )
        for topo_name in ["bond","angle","dihedral","improper"]:
            self_topo = self._get_topo_list(topo_name)
            second_topo = world_2._get_topo_list(topo_name)
            second_topo = [Topology(*[a+atom_offset for a in t.get_atoms()]  # Each atom ID referenced in new topo is offset
                ,type_id=t.type_id) for t in second_topo]
            self_topo.extend(second_topo)  # Mutability lets us extend the original list.


    def add_atom(self,type_id,charge,pos,img_flags = [0,0,0]):
        """
        Add atom to system
        
        Parameters
        ----------
        type_id : int
            The atom type index
        charge : float
            Atom partial charge
        pos : list of float
            List of length 3, defining x, y and z position of atom respectively
        img_flags : list of int, optional
            List of length 3, defining x, y and z image flags of atom respectively (default is [0,0,0])
        """
        # TODO: Validate pos & image flags?
        # Or convert to absolute pos, back to image flags on output?

        if not 0 <= int(type_id) < self.n_atom_types:
            raise ValueError(f"Atom type {type_id} is out of bounds")

        self.atoms.append(Atom(type_id,charge,pos,img_flags))
 
    def _add_topo(self, topo_name = None, *atoms, type_id=0):
        topo_list = self._get_topo_list(topo_name)

        if type_id >= len(self._get_topo_param_list(topo_name)):
            raise ValueError(f"{topo_name.capitalize()} type {type_id} is out of bounds.")

        # verify valid topology data
        new_topo = Topology(*self._validate_topo_atoms(*atoms), type_id=type_id)
        
        # Check duplicate
        if self._get_topo_index(topo_name, *new_topo) >= 0:
            print(f"WARNING: {topo_name} {new_topo} already exists.")
            return

        topo_list.append(new_topo)

    def add_bond(self, atom1, atom2, type_id=0):
        self._add_topo("bond", atom1, atom2, type_id=type_id)
    
    def add_angle(self, atom1, atom2, atom3, type_id=0):
        self._add_topo("angle", atom1, atom2, atom3, type_id=type_id)

    def add_dihedral(self, atom1, atom2, atom3, atom4, type_id=0):
        self._add_topo("dihedral", atom1, atom2, atom3, atom4, type_id=type_id)

    def add_improper(self, atom1, atom2, atom3, atom4, type_id=0):
        self._add_topo("improper", atom1, atom2, atom3, atom4, type_id=type_id)

    # Delete atom and all associated topology (bonds, angles, ...)
    def delete_atom(self,atom_index):
        """
        Delete atom. Automatically removes all topology referencing the deleted atom.

        Parameters
        ----------
        atom_index : int
            Index of atom to be deleted.

        Returns
        -------
        Atom
            Atom object that was removed from the system.
        """
        
        # Find all topology containing atom to be deleted and remove it
        topo_deleted = False
        for topo_name in ["bond", "angle", "dihedral", "improper"]:
            topo_indices = self._find_connected_topo_indices(topo_name, atom_index)
            if len(topo_indices) > 0: topo_deleted = True
            self._delete_topos(topo_name, *topo_indices)
        
        if topo_deleted:
            print(f"WARNING: Deleting atom {atom_index} has led to the deletion of topology.")

        # Any atoms with higher ID than deleted atom will have their ID changed
        # So we must update all topology to reference the correct atoms
        self._shift_topo_indices(range(atom_index+1,self.count_atoms()), shift = -1)

        return self.atoms.pop(atom_index)

    
    def _delete_topos(self, topo_name, *delete_indices):
        topo_list = self._get_topo_list(topo_name)
        new_topo = [topo for i,topo in enumerate(topo_list) if i not in delete_indices]
        self._set_topo_list(topo_name, new_topo)

    def delete_bonds(self, *delete_indices):
        self._delete_topos("bond", *delete_indices)

    def delete_angles(self, *delete_indices):
        self._delete_topos("angle", *delete_indices)
        
    def delete_dihedrals(self,*delete_indices):
        self._delete_topos("dihedral", *delete_indices)

    def delete_impropers(self, *delete_indices):
        self._delete_topos("improper", *delete_indices)

    def get_isolated_atoms(self, atom_indices : List[int], return_equivalences : bool = False) -> 'World':
        """
        Return a new world containing only the atoms in atom_indices, and relevant topology.
        
        Parameters
        ----------
        atom_indices : List[int]
            The indices of the atoms to be included in the returned World.
        return_equivalences : bool
            If true, a dictionary containing the conversions from all atoms to fragment atoms IDs is returned also.

        Returns
        -------
        World
            A world containing only the specified atoms and relevant topology.
        """

        try:
            atom_indices = sorted([int(a) for a in atom_indices])  # Sort so things work properly later
        except:
            raise TypeError(f"atom_indices must be an Iterable object containing integers, not {type(atom_indices)}.")

        equivalences = dict()

        # Create new world object with the same dimensions as the world the atoms are being taken from.
        isolated = World(self.xlo, self.ylo, self.zlo, self.xhi, self.yhi, self.zhi)

        isolated.atom_type_params = self.atom_type_params
        for topo_name in ("bond","angle","dihedral","improper"):
            isolated_topo_params = isolated._get_topo_param_list(topo_name)
            # Since isolated was only just made, isolated_topo_params will be empty. Therefore, .extend() sets the whole list.
            isolated_topo_params.extend(self._get_topo_param_list(topo_name))

        # Copy all atoms of interest as independent objects.
        isolated.atoms = [deepcopy(a) for i,a in enumerate(self.atoms) if i in atom_indices]
        # Kept atoms will be kept in the same order, but have sequential indices - map lets us convert topo references.
        equivalences.update({a:i for i,a in enumerate(atom_indices)})  # Update empty dictionary argument will update dictionary outside function.
        for topo_name in ("bond","angle","dihedral","improper"):
            isolated_topo = isolated._get_topo_list(topo_name)
            self_topo = self._get_topo_list(topo_name)

            # First line  - Atom IDs change when atoms are isolated. Use atom_id_map to make topology referencing new IDs.
            # Second line - Only copy topology items if *all* of the referenced atoms are being isolated.
            relevant_topo = [Topology(*[equivalences[a] for a in t], type_id = t.type_id)
                                for t in self_topo if not any(a not in atom_indices for a in t)]
            isolated_topo.extend(relevant_topo)

        if return_equivalences: return isolated, equivalences
        else: return isolated


    def _get_connected_topo(self, topo_name : str, atom_index : int) -> List:
        """Returns list of topology objects of given type that include atom_index."""
        topo_list = self._get_topo_list(topo_name)
        return [t for t in topo_list if atom_index in t]

    def get_connected_bonds(self, atom_index : int) -> List:
        """Returns list of bond objects that include atom_index."""
        return self._get_connected_topo("bond",atom_index)

    def get_connected_angles(self, atom_index : int) -> List:
        """Returns list of angle objects that include atom_index."""
        return self._get_connected_topo("angle",atom_index)

    def get_connected_dihedrals(self, atom_index : int) -> List:
        """Returns list of dihedral objects that include atom_index."""
        return self._get_connected_topo("dihedral",atom_index)

    def get_connected_impropers(self, atom_index : int) -> List:
        """Returns list of improper objects that include atom_index."""
        return self._get_connected_topo("improper",atom_index)
    
    def _find_connected_topo_indices(self, topo_name : str, atom_index : int) -> List:
        """Returns indices of topology objects containing atom"""
        topo_list = self._get_topo_list(topo_name)
        return [i for i,t in enumerate(topo_list) if atom_index in t]

    def find_connected_bond_indices(self, atom_index):
        """Returns indices of bonds containing atom"""
        return self._find_connected_topo_indices("bond", atom_index)

    def find_connected_angle_indices(self, atom_index) -> List[int]:
        """Returns indices of angles containing atom"""
        return self._find_connected_topo_indices("angle", atom_index)

    def find_connected_dihedral_indices(self, atom_index):
        """Returns indices of dihedrals containing atom"""
        return self._find_connected_topo_indices("dihedral", atom_index)

    def find_connected_improper_indices(self, atom_index):
        """Returns indices of impropers containing atom"""
        return self._find_connected_topo_indices("improper", atom_index)
    
    def _find_topo_connected_atoms(self, topo_name : str, atom_index:int):
        """
        Get atom IDs of all atoms connected to reference atom 'atom_index' by topology 'topo_name'.
        
        Arguments
        ---------
        topo_name : str
            String-type name of topology of interest. Allowed values are 'bond', 'angle', 'dihedral' and 'improper'.
        atom_index : int
            Atom ID number of reference atom
        """

        topo = self._get_connected_topo(topo_name,atom_index)

        # Conversion to a set and back removes duplicates - feels bad.
        return list(set([a for t in topo for a in t if a != atom_index]))

    def find_bonded_atoms(self,atom_index : int) -> List[int]:
        return self._find_topo_connected_atoms("bond", atom_index)
    
    def find_angle_atoms(self,atom_index : int) -> List[int]:
        return self._find_topo_connected_atoms("angle", atom_index)

    def find_dihedral_atoms(self,atom_index : int) -> List[int]:
        return self._find_topo_connected_atoms("dihedral", atom_index)

    def find_improper_atoms(self,atom_index : int) -> List[int]:
        return self._find_topo_connected_atoms("improper", atom_index)

    def _get_topo_index(self, topo_name, *atoms):
        """
        Check if topology exists between atom1 and atom2.
        returns -1 if doesn't exist, topo index if does exist.

        BE CAREFUL using output as list index: -1 is a valid index in python
        """

        # Get list of relevant topology
        topo_list = self._get_topo_list(topo_name)

        # Verify atoms make valid topology and are ordered correctly
        atoms = self._validate_topo_atoms(*atoms)

        for i,t in enumerate(topo_list):
            if atoms == t:
                return i  # Match found, return existing topo index
        return -1  # If we never found a match, topo doesn't exist, return -1

    def get_bond_index(self, atom1, atom2):
        return self._get_topo_index("bond", atom1, atom2)
    
    def get_angle_index(self, atom1, atom2, atom3):
        return self._get_topo_index("angle", atom1, atom2, atom3)

    def get_dihedral_index(self, atom1, atom2, atom3, atom4):
        return self._get_topo_index("dihedral", atom1, atom2, atom3, atom4)

    def get_improper_index(self, atom1, atom2, atom3, atom4):
        return self._get_topo_index("improper", atom1, atom2, atom3, atom4)

    def count_atoms(self):
        return len(self.atoms)

    def append_atom_types(self, new_atom_types : list) -> None:
        try:
            if len(new_atom_types) == 0:
                raise ValueError("List of new atom types must not be empty!")
        except:
            raise TypeError(f"new_atom_types must be a list or other iterable, not {type(new_atom_types)}.")

        try:
            new_atom_types = [float(m) for m in new_atom_types]
        except ValueError as err:
            raise Exception(f"One or more atom type properties is invalid in {new_atom_types}.") from err

        self.atom_type_params.extend(new_atom_types)


    def _append_topo_types(self, topo_name : str, new_topo_types : list) -> None:
        try:
            if len(new_topo_types) == 0:
                raise ValueError("List of new topology types must not be empty!")
        except:
            raise TypeError(f"new_topo_types must be a list or other iterable, not {type(new_topo_types)}.")
        
        param_list = self._get_topo_param_list(topo_name)

        param_list.extend(new_topo_types)

    def append_bond_types(self, new_bond_types : list) -> None:
        self._append_topo_types("bond", new_bond_types)
    
    def append_angle_types(self, new_angle_types : list) -> None:
        self._append_topo_types("angle", new_angle_types)
    
    def append_dihedral_types(self, new_dihedral_types : list) -> None:
        self._append_topo_types("dihedral", new_dihedral_types)
    
    def append_improper_types(self, new_improper_types : list) -> None:
        self._append_topo_types("improper", new_improper_types)
    
        
    def set_atom_type_mass(self, type_index : int, mass : float) -> None:
        try:
            type_index = int(type_index)
        except ValueError as err:
            raise Exception(f"Invalid type for atom type index : {type(type_index)}")

        try:
            mass = float(mass)
        except ValueError as err:
            raise Exception(f"Invalid type for atom type mass : {type(mass)}")

        if mass < 0:
            raise ValueError(f"Atom type mass must be positive or zero, not {mass}.")

        self.atom_type_params[type_index] = mass

    def add_atom_type(self, mass, type_id = "append"):
        """Add or insert a new atom type. If type_id is not specifies, appends the new value."""

        try:
            mass = float(mass)
        except:
            raise TypeError(f"Atom mass must be a float, not {type(mass)}.")

        try:
            if type_id != "append": type_id = int(type_id)
        except:
            raise TypeError(f"Atom type ID must be an integer or None, not {type(type_id)}.")

        if type_id == "append" or type_id == len(self.atom_type_params):
            self.atom_type_params.append(mass)
        elif type_id < len(self.atom_type_params):
            self.atom_type_params.insert(type_id,mass)
            self._shift_atom_types(range(type_id,len(self.atom_type_params)),shift=+1)
        else:
            raise ValueError(f"type_id must be between 0 and the number of existing atom types ({self.n_atom_types}), inclusive. Received {type_id}.")

    def _add_topo_type(self, topo_name, *params, topo_id = "append"):
        """
        Add a bond, angle, ... type, with given interaction parameters.

        If topo_id is unspecified, appends new type.

        If topo_id is an integer, inserts topo type at given index.
        """

        # Get list of relevant topology types
        param_list = self._get_topo_param_list(topo_name)

        # Check topo_id has valid value
        if topo_id == "append":
            topo_id = len(param_list)
        elif topo_id not in range(0, len(param_list)+1):
            raise ValueError(f"{topo_name.capitalize()} type {topo_id} is out of bounds.")
        
        # Offset topology type references in topology list
        self._shift_topo_types(topo_name, range(topo_id,len(param_list)), shift = +1)

        # Insert bond type at relevant index
        param_list.insert(topo_id, params)
        

    def add_bond_type(self, *params, bond_id = "append"):
        self._add_topo_type("bond", *params, topo_id = bond_id)
    
    def add_angle_type(self, *params, angle_id = "append"):
        self._add_topo_type("angle", *params, topo_id = angle_id)

    def add_dihedral_type(self, *params, dihedral_id = "append"):
        self._add_topo_type("dihedral", *params, topo_id = dihedral_id)
    
    def add_improper_type(self, *params, improper_id = "append"):
        self._add_topo_type("improper", *params, topo_id = improper_id)

    
    def _validate_topo_atoms(self, *atoms):
        """Used to check bond/angle/dihedral atom inputs are valid and order them properly"""
        # Check atom indices are integers
        try:
            atoms = [int(a) for a in atoms]
        except ValueError as err:
            raise Exception(f"One or more atom indices in {atoms} is not an integer.") from err

        # Ensure first atom index is smaller than last
        if atoms[0] > atoms[-1]:
            atoms = atoms[::-1]  # Reverses tuple of atoms

        for at in atoms:
            if at >= self.count_atoms():
                raise IndexError(f"Atom index {at} is out of bounds.")

        return atoms

    def _shift_atom_types(self, target_type_range, shift : int) -> None:
        """Shift the atom type numbers of atoms in the world. Useful for when atom types are added/removed within the list.
        
        Parameters
        ----------
        target_type_range : Iterable (list, range, ...)
            Pre-shift atom type indices that need to be shifted.
        shift : int
            Amount to shift selected indices by.

        Returns
        -------
        None
        """

        if not isinstance(shift, int):
            raise TypeError(f"Shift must be of type int, not {type(shift)}.")
        try:
            target_type_range = [int(t) for t in target_type_range]
        except:
            raise TypeError(f"target_type_range must be a list of integers, not {target_type_range}.")
        for a in self.atoms:
            if a.type_id in target_type_range:
                a.type_id += shift

    def _shift_topo_atom_indices(self, target_atom_range, shift : int) -> None:
        """
        Shift the atom indices referenced in topology. Useful for when atoms are removed/added anywhere but the end of the list.

        Parameters
        ----------
        target_atom_range : Iterable (list, range, ...)
            Pre-shift atom indices that need to be shifted.
        shift : int
            Amount to shift selected indices by.

        Returns
        -------
        None
        """

        if not isinstance(shift, int):
            raise TypeError(f"Shift must be of type int, not {type(shift)}.")
        for topo_list in [self.bonds, self.angles, self.dihedrals, self.impropers]:
            for topo in topo_list:
                for i,atom in enumerate(topo):
                    if atom in target_atom_range: topo[i] += shift
                

    def _shift_topo_types(self, topo_name, target_type_range : Iterable, shift : int):
        """
        Shift the topology type indices referenced in topology.

        Useful for when topology types are removed/added anywhere but the end of the list.

        Parameters
        ----------
        target_type_range : Iterable
            Pre-shift type indices that need to be shifted.
        shift : int
            Amount to shift type indices by.
        """

        if not isinstance(shift, int):
            raise TypeError(f"Shift must be of type int, not {type(shift)}.")
        for topo in self._get_topo_list(topo_name):
            if topo.type_id in target_type_range:
                topo.type_id += shift
            # TODO: Check for out-of-bounds types

    # Return list of all the topology of name topo_name, where topo_name is a string
    # Allowed values are: bond, angle, dihedral, or improper
    def _get_topo_list(self, topo_name):
        """
        Return list of all the topology of name topo_name.

        Parameters
        ----------
        topo_name : str
            Name of topology type. Allowed values are: bond, angle, dihedral, or improper.
        """

        if not isinstance(topo_name,str):
            raise TypeError(f"Topo name must be a string, not {type(topo_name)}.")

        if topo_name == "bond": return self.bonds
        elif topo_name == "angle": return self.angles
        elif topo_name == "dihedral": return self.dihedrals
        elif topo_name == "improper": return self.impropers
        else: raise ValueError(f"{repr(topo_name)} is not a valid value of topo_name.")
    
    # Set list of topology to supplied list topo_list
    def _set_topo_list(self, topo_name : str, topo_list : list):
        if not isinstance(topo_name,str):
            raise TypeError(f"Topo name must be a string, not {type(topo_name)}.")

        if topo_name == "bond": self.bonds = topo_list
        elif topo_name == "angle": self.angles = topo_list
        elif topo_name == "dihedral": self.dihedrals = topo_list
        elif topo_name == "improper": self.impropers = topo_list
        else: raise ValueError(f"{repr(topo_name)} is not a valid value of topo_name.")

    def _get_topo_param_list(self, topo_name):
        if not isinstance(topo_name,str):
            raise TypeError(f"Topo name must be a string, not {type(topo_name)}.")

        if topo_name == "bond": return self.bond_params
        elif topo_name == "angle": return self.angle_params
        elif topo_name == "dihedral": return self.dihedral_params
        elif topo_name == "improper": return self.improper_params
        else: raise ValueError(f"{repr(topo_name)} is not a valid value of topo_name.")
 