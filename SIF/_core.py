"""
---Simulation Input Finagler Core---
Defining classes and methods for
    loading, manipulating and saving
    simulation info.
"""

from argparse import ArgumentError
from typing import Iterable,List,Union
import math
from copy import deepcopy
import logging
import numpy as np
import itertools

# Logger object, for nicely formatting warnings, errors etc.
logger = logging.getLogger(__name__)

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
    pos : array of float
        Numpy array of shape (3,), defining x, y and z position of atom respectively
    img_flags : array of int
        Numpy array of shape (3,), defining x, y and z image flags of atom respectively
    """

    # Initialiser, used when making new atom objects.
    def __init__(self, type_id:int, charge:float, pos:np.ndarray, img_flags:np.ndarray=np.zeros(3)):

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
            pos = np.array(pos)
        except ValueError as err:
            raise Exception(f"Atom position {pos} should be something that can be turned into a numpy array.") from err
        
        # Check pos has 3 elements
        try:
            assert pos.shape == (3,)
        except AssertionError as err:
            raise Exception(f"Atom position {pos} should have exactly 3 elements.") from err

        # Check pos is all floats
        try:
            pos = pos.astype(np.float64)
        except ValueError as err:
            raise Exception(f"Atom position {pos} contains a non-float element.") from err

        # Check img_flags is a list-type variable
        try:
            img_flags = np.array(img_flags)
        except ValueError as err:
            raise Exception(f"Atom image flags {img_flags} should be something that can be turned into a numpy array.") from err

        # Check img_flags has 3 elements
        try:
            assert img_flags.shape == (3,)
        except AssertionError as err:
            raise Exception(f"Atom image flags {img_flags} should have exactly 3 elements.") from err

        # Check img_flags is all ints
        try:
            img_flags = img_flags.astype(np.int32)
        except ValueError as err:
            raise Exception(f"Atom image flags {img_flags} contains a non-integer element.") from err


        # Attributes
        self.type_id = type_id
        self.charge = charge
        self.pos = pos
        self.img_flags = img_flags
    
    def __copy__(self):
        """Ensures pos and image flag lists are copied properly too."""
        return Atom(self.type_id,self.charge,np.copy(self.pos),np.copy(self.img_flags))
    
    def copy(self):
        """Return an identical copy of the atom"""
        return self.__copy__()


class AtomType:
    """
    Information on a type of atom defined in a system.
    
    Attributes
    ----------
    mass : float
        Mass of the given type of atom
    name : string
        Unique atom type name (often defined by force fields)
    """

    def __init__(self, mass:float, name:str=None):
        self.mass=float(mass)
        self.name=str(name)


class Topology(list):
    """
    A class used to represent a topology instance (a bond, an angle, ...), inheriting from list.

    By inheriting list, we get to say Topology[i] to get atom i of the topology.

    Attributes
    ----------
    type_id : int
        The type index of the given topology (used to define interaction parameters)
    name : string
        Unique type name indicating the type of topology described.
    
    Methods
    -------
    get_atoms():
        Returns an ordered list of atoms defining the topology
    """
    def __init__(self, *atoms, type_id:int=0):
        self.type_id = int(type_id)
        super().__init__(atoms)  # Initialise the list part
    
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

class TopologyType:
    """Similar to AtomType, stores the parameters and name of a type of topology interaction.
    
    Attributes
    ----------
    parameters : list[float]
        Interaction parameters defined for given topology type. Forcefield-specific.
    name : string (Default = None)
        Topology type name used to uniquely and consistently identify the interaction.
    """

    def __init__(self,*params,name:str=None):
        if not isinstance(params,Iterable):
            params = [params]
        try:
            params = [float(p) for p in params]
        except:
            raise TypeError("Some topology interaction parameters were not floats!")

        self.parameters=params
        self.name = str(name)

    def __copy__(self:'TopologyType'):
        return TopologyType(*self.parameters,name=self.name)
    
    def copy(self:'TopologyType'):
        """Return an identical but separate copy of this object."""
        return self.__copy__()

class World:

    # List of directly implemented topology type names
    #   using other types will probably start breaking things!
    # Not sure how to get around this but shouldn't usually be
    #   a problem.
    _available_topo_types = ("bond","angle","dihedral","improper")
    
    def __init__(self, xlo=0., ylo=0., zlo=0., xhi=0., yhi=0., zhi=0.):
        
        # Box Dimensions
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi

        # Atoms & topology
        self.atoms : List[Atom] = []
        self.bonds : List[Topology] = []  
        self.angles : List[Topology] = []
        self.dihedrals : List[Topology] = []
        self.impropers : List[Topology] = []

        # Force field data
        self.atom_types : List[AtomType] = []  # atom_types contains list of AtomType objects
        # Topology types contain lists of TopologyType objects
        self.bond_types : List[TopologyType] = []
        self.angle_types : List[TopologyType] = []
        self.dihedral_types : List[TopologyType] = []
        self.improper_types : List[TopologyType] = []

        # Additional metadata
        self.fragments = {}  # Dict[str,list[int]]. Used by LAMMPS molecule files.

    @property
    def n_atom_types(self):
        """Return the number of atom types available in the system."""
        return len(self.atom_types)
    
    @property
    def n_atoms(self):
        """The number of atoms in the system."""
        return len(self.atoms)

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

    def add_fragment(self:'World', frag_id : str = "", atom_ids : List[int] = []):
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
        frag_id = str(frag_id)

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
    
    def set_dims(self:'World', xlo=None, ylo=None, zlo=None, xhi=None, yhi=None, zhi=None):
        """Set World bounding box dimensions."""
        if xlo: self.xlo = float(xlo)
        if xhi: self.xhi = float(xhi)
        if ylo: self.ylo = float(ylo)
        if yhi: self.yhi = float(yhi)
        if zlo: self.zlo = float(zlo)
        if zhi: self.zhi = float(zhi)

    def merge(self:'World', world_2:'World', type_merge_style="Overwrite") -> None:
        """
        Merge a second world into this world, arranging topology and atoms automatically.
        Returns nothing, instead directly modifying this world.

        Arguments
        ---------
        world_2 : World
            The world to merge into this one.
        type_merge_style : string
            How to merge the atom and topology type information of the two worlds.
            - "Keep": Do not overwrite any types, keep those from this world.
            - "Overwrite" - Naive method, overwrites the atom types if world_2 has more types of atoms. (Default)
            - "Merge" - Uses atom type names to merge groups without duplicates.
            - "Append" - ***NOT IMPLEMENTED***
        """
        # Sort atom/topo types first
        atom_offset = self.n_atoms
        n_params = [{   "atoms": self.n_atom_types,
                        "bonds": len(self.bond_types),
                        "angles": len(self.angle_types),
                        "dihedrals": len(self.dihedral_types),
                        "impropers": len(self.improper_types)},
                    {   "atoms": world_2.n_atom_types,
                        "bonds": len(world_2.bond_types),
                        "angles": len(world_2.angle_types),
                        "dihedrals": len(world_2.dihedral_types),
                        "impropers": len(world_2.improper_types)}]

        unequal_types = [n_params[0][t] != n_params[1][t] for t in ("atoms", "bonds", "angles", "dihedrals", "impropers")]
        if True in unequal_types: logger.warning("Merged systems have different numbers of atom and/or topology types.")
        w2_atoms = deepcopy(world_2.atoms)

        if type_merge_style == "Overwrite":
            if n_params[1]["atoms"] > n_params[0]["atoms"]: self.atom_types = world_2.atom_types.copy()
            if n_params[1]["bonds"] > n_params[0]["bonds"]: self.bond_types = world_2.bond_types.copy()
            if n_params[1]["angles"] > n_params[0]["angles"]: self.angle_types = world_2.angle_types.copy()
            if n_params[1]["dihedrals"] > n_params[0]["dihedrals"]: self.dihedral_types = world_2.dihedral_types.copy()
            if n_params[1]["impropers"] > n_params[0]["impropers"]: self.improper_types = world_2.improper_types.copy()
        elif type_merge_style == "Merge":
            new_types = []
            type_names_w1 = [t.name for t in self.atom_types]
            type_names_w2 = [t.name for t in world_2.atom_types]
            if None in type_names_w1 or None in type_names_w2:
                raise ValueError("Merge style can't be used if not every atom has a type name!")
            w2_i = 0
            for type_w1 in self.atom_types:
                type_w2 = world_2.atom_types[w2_i]
                if type_w1.name not in type_names_w2:
                    new_types.append(type_w1)
                elif type_w1.name == type_w2.name:
                    if type_w1.mass != type_w2.mass:
                        raise ValueError("Atom types with the same name but different parameters are present in merging worlds!")
                    new_types.append(type_w1)
                    w2_i += 1
                else:  # Take all w2 types up to and including the one matching w1_i.
                    w2_match_i = type_names_w2.index(type_w1.name)
                    new_types.extend(world_2.atom_types[w2_i:w2_match_i+1])
                    w2_i = w2_match_i + 1
            new_types.extend(world_2.atom_types[w2_i:])

            # Manually looking up each one is quite bad... Could have built a lookup?
            for atom in self.atoms:
                type_name = self.atom_types[atom.type_id]
                atom.type_id = new_types.index(type_name)
            for atom in w2_atoms:
                type_name = world_2.atom_types[atom.type_id]
                atom.type_id = new_types.index(type_name)
                    
        self.atoms.extend(
            [a.copy() for a in world_2.atoms]
            )
        for topo_name in World._available_topo_types:
            self_topo = self._get_topo_list(topo_name)
            second_topo = world_2._get_topo_list(topo_name)
            second_topo = [Topology(*[a+atom_offset for a in t.get_atoms()]  # Each atom ID referenced in new topo is offset
                ,type_id=t.type_id) for t in second_topo]
            self_topo.extend(second_topo)  # Because we have a *reference* to the list, this will extend the original.

    def _topo_type_id_from_name(self:'World',topo_name:str,type_name:str)->int:
        """Return the (first, and hopefully only) index (type_id) of the topology type with name type_name."""
        topo_types = self._get_topo_type_list(topo_name)
        return next(i for i,t in enumerate(topo_types) if t.name == type_name)
    
    def bond_type_id_from_name(self:'World', type_name:str)->int:
        """Return the (first, and hopefully only) index (type_id) of the bond type with name type_name."""
        return self._topo_type_id_from_name("bond",type_name)
    
    def angle_type_id_from_name(self:'World', type_name:str)->int:
        """Return the (first, and hopefully only) index (type_id) of the angle type with name type_name."""
        return self._topo_type_id_from_name("angle",type_name)
    
    def dihedral_type_id_from_name(self:'World', type_name:str)->int:
        """Return the (first, and hopefully only) index (type_id) of the dihedral type with name type_name."""
        return self._topo_type_id_from_name("dihedral",type_name)
    
    def improper_type_id_from_name(self:'World', type_name:str)->int:
        """Return the (first, and hopefully only) index (type_id) of the improper type with name type_name."""
        return self._topo_type_id_from_name("improper",type_name)
    
    def atom_type_id_from_name(self:'World', type_name:str)->int:
        """Return the (first, and hopefully only) index (type_id) of the atom type with name type_name."""
        return next(i for i,t in enumerate(self.atom_types) if t.name == type_name)


    def add_atom(self:'World',
                 type_id:Union[int,str]=None,
                 charge:float=0,
                 pos:np.ndarray=np.zeros(3),
                 img_flags:np.ndarray=np.zeros(3)):
        """
        Add an atom to the system.
        
        Parameters
        ----------
        type_id : int or str
            The atom type index (int) or name (str)
        charge : float
            Atom partial charge (default is 0)
        pos : list of float
            List of length 3, defining x, y and z position of atom respectively (default is [0,0,0])
        img_flags : list of int, optional
            List of length 3, defining x, y and z image flags of atom respectively (default is [0,0,0])
        """
        # TODO: Validate pos & image flags?
        # Or convert to absolute pos, back to image flags on output?

        if isinstance(type_id, str):
            type_id = self.atom_type_id_from_name(type_id)
        elif not 0 <= int(type_id) < self.n_atom_types:
            raise ValueError(f"Atom type {type_id} is out of bounds")

        self.atoms.append(Atom(type_id,charge,pos,img_flags))
 
    def _add_topo(self:'World', topo_name, *atoms, type_id:Union[int,str]=None):
        """Add a new topology of any type between any number of atoms (specified by index) to the system."""

        assert(type_id is not None, "type_id is a required argument!")

        if isinstance(type_id, str):
            type_id = self._topo_type_id_from_name(topo_name, type_id)

        topo_list = self._get_topo_list(topo_name)

        if type_id >= len(self._get_topo_type_list(topo_name)):
            raise ValueError(f"{topo_name.capitalize()} type {type_id} is out of bounds.")

        # verify valid topology data
        new_topo = Topology(*self._validate_topo_atoms(*atoms), type_id=type_id)
        
        # Check duplicate
        if self._get_topo_index(topo_name, *new_topo) is not None:
            logger.warning(f"{topo_name} {new_topo} already exists.")
            return

        topo_list.append(new_topo)

    def add_bond(self, atom1:int, atom2:int,  type_id:Union[int,str]=None):
        """Create a bond interaction between atom indices atom1 and atom2."""
        self._add_topo("bond", atom1, atom2, type_id=type_id)
    
    def add_angle(self, atom1:int, atom2:int, atom3:int, type_id:Union[int,str]=None):
        """Create an angle interaction between atom indices atom1, atom2 and atom3."""
        self._add_topo("angle", atom1, atom2, atom3, type_id=type_id)

    def add_dihedral(self, atom1:int, atom2:int, atom3:int, atom4:int, type_id:Union[int,str]=None):
        """Create a dihedral interaction between atom indices atom1, atom2, atom3 and atom4."""
        self._add_topo("dihedral", atom1, atom2, atom3, atom4, type_id=type_id)

    def add_improper(self, atom1:int, atom2:int, atom3:int, atom4:int, type_id:Union[int,str]=None):
        """Create an improper dihedral interaction between atom indices atom1, atom2, atom3 and atom4."""
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
        for topo_name in World._available_topo_types:
            topo_indices = self._find_connected_topo_indices(topo_name, atom_index)
            if len(topo_indices) > 0: topo_deleted = True
            self._delete_topos(topo_name, *topo_indices)
        
        if topo_deleted:
            logger.warning(f"Deleting atom {atom_index} has led to the deletion of topology.")

        # Any atoms with higher ID than deleted atom will have their ID changed
        # So we must update all topology to reference the correct atoms
        self._shift_topo_atom_indices(range(atom_index+1,self.n_atoms), shift = -1)

        return self.atoms.pop(atom_index)

    
    def _delete_topos(self, topo_name, *delete_indices):
        """Delete topology of type topo_name specified by index. Shifts topology indices to keep them sequential."""
        topo_list = self._get_topo_list(topo_name)
        new_topo = [topo for i,topo in enumerate(topo_list) if i not in delete_indices]
        self._set_topo_list(topo_name, new_topo)

    def delete_bonds(self, *delete_indices):
        """Delete bonds specified by index. Shifts bond indices to keep them sequential."""
        self._delete_topos("bond", *delete_indices)

    def delete_angles(self, *delete_indices):
        """Delete angles specified by index. Shifts angle indices to keep them sequential."""
        self._delete_topos("angle", *delete_indices)
        
    def delete_dihedrals(self,*delete_indices):
        """Delete dihedrals specified by index. Shifts dihedral indices to keep them sequential."""
        self._delete_topos("dihedral", *delete_indices)

    def delete_impropers(self, *delete_indices):
        """Delete improper dihedrals specified by index. Shifts improper dihedral indices to keep them sequential."""
        self._delete_topos("improper", *delete_indices)

    def get_isolated_atoms(self, atom_indices : List[int], return_equivalences : bool = False) -> 'World':
        """
        Return a new world containing only the atoms in atom_indices, and relevant topology.
        
        Parameters
        ----------
        atom_indices : List[int]
            The indices of the atoms to be included in the returned World.
        return_equivalences : bool (Default = False)
            If true, a dictionary containing the conversions from origin world atom IDs to isolated atom IDs is returned also.

        Returns
        -------
        World
            A world containing only the specified atoms and relevant topology.
        Dict[int,int] (optional - only if return_equivalences = True)
            A dictionary to convert from original world index (key) to isolated world index (value).
        """

        try:
            atom_indices = sorted([int(a) for a in atom_indices])  # Sort so things work properly later
        except:
            raise TypeError(f"atom_indices must be an Iterable object containing integers, not {type(atom_indices)}.")

        equivalences = dict()

        # Create new world object with the same dimensions as the world the atoms are being taken from.
        isolated = World(self.xlo, self.ylo, self.zlo, self.xhi, self.yhi, self.zhi)

        isolated.atom_types = deepcopy(self.atom_types)
        for topo_name in World._available_topo_types:
            isolated_topo_params = isolated._get_topo_type_list(topo_name)
            # Since isolated was only just made, isolated_topo_params will be empty. Therefore, .extend() sets the whole list.
            isolated_topo_params.extend(self._get_topo_type_list(topo_name))

        # Copy all atoms of interest as independent objects.
        isolated.atoms = [deepcopy(a) for i,a in enumerate(self.atoms) if i in atom_indices]
        # Kept atoms will be kept in the same order, but have sequential indices - map lets us convert topo references.
        equivalences.update({a:i for i,a in enumerate(atom_indices)})  # Update empty dictionary argument will update dictionary outside function.
        for topo_name in World._available_topo_types:
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
        Get atom IDs of all atoms directly connected to reference atom 'atom_index' by topology 'topo_name'.
        
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
        """
        Get atom IDs of all atoms directly connected to reference atom 'atom_index' by bonds.
        
        Arguments
        ---------
        atom_index : int
            Atom ID number of reference atom
        """
        return self._find_topo_connected_atoms("bond", atom_index)
    
    def find_angle_atoms(self,atom_index : int) -> List[int]:
        """
        Get atom IDs of all atoms directly connected to reference atom 'atom_index' by angle topologies.
        
        Arguments
        ---------
        atom_index : int
            Atom ID number of reference atom
        """
        return self._find_topo_connected_atoms("angle", atom_index)

    def find_dihedral_atoms(self,atom_index : int) -> List[int]:
        """
        Get atom IDs of all atoms directly connected to reference atom 'atom_index' by dihedral topologies.
        
        Arguments
        ---------
        atom_index : int
            Atom ID number of reference atom
        """
        return self._find_topo_connected_atoms("dihedral", atom_index)

    def find_improper_atoms(self,atom_index : int) -> List[int]:
        """
        Get atom IDs of all atoms directly connected to reference atom 'atom_index' by improper dihedral topologies.
        
        Arguments
        ---------
        atom_index : int
            Atom ID number of reference atom
        """
        return self._find_topo_connected_atoms("improper", atom_index)

    def _get_topo_index(self, topo_name, *atoms):
        """
        Check if topology object exists linking given atom indices together.
        Returns None if doesn't exist, topo index if does exist.
        """

        # Get list of relevant topology
        topo_list = self._get_topo_list(topo_name)

        # Verify atoms make valid topology and are ordered correctly
        atoms = self._validate_topo_atoms(*atoms)

        for i,t in enumerate(topo_list):
            if atoms == t:
                return i  # Match found, return existing topo index
        return None  # If we never found a match, topo doesn't exist, return -1

    def get_bond_index(self, atom1, atom2):
        """
        Check if a bond exists between atom indices atom1 and atom2.
        Returns None if doesn't exist, bond index if does exist.
        """
        return self._get_topo_index("bond", atom1, atom2)
    
    def get_angle_index(self, atom1, atom2, atom3):
        """
        Check if an angle topology exists between atom indices atom1, atom2 and atom3.
        Returns None if doesn't exist, angle index if does exist.
        """
        return self._get_topo_index("angle", atom1, atom2, atom3)

    def get_dihedral_index(self, atom1, atom2, atom3, atom4):
        """
        Check if a dihedral topology exists between atom indices atom1, atom2, atom3 and atom4.
        Returns None if doesn't exist, dihedral index if does exist.
        """
        return self._get_topo_index("dihedral", atom1, atom2, atom3, atom4)

    def get_improper_index(self, atom1, atom2, atom3, atom4):
        """
        Check if an improper dihedral topology exists between atom indices atom1, atom2, atom3 and atom4.
        Returns None if doesn't exist, improper dihedral index if does exist.
        """
        return self._get_topo_index("improper", atom1, atom2, atom3, atom4)

    def count_atoms(self):
        """DEPRECATED: Use the n_atoms property."""
        logger.warning("count_atoms is deprecated. Use the n_atoms property.")
        return len(self.n_atoms)

    def append_atom_types(self, new_atom_types : List[AtomType]) -> None:
        """Define new atom types, appended to the end of the current list.
        
        Arguments
        ---------
        new_atom_types : List[AtomType]
            List of new AtomTypes to be added to the system.
        """
        try:
            if len(new_atom_types) == 0:
                logger.warning("Just appended a new atom type list of length 0!")
                return
        except:
            raise TypeError(f"new_atom_types must be a list or other iterable, not {type(new_atom_types)}.")

        # If the list contains anything that isn't an AtomType object, throw an error.
        if False in [isinstance(t,AtomType) for t in new_atom_types]:
            raise TypeError(f"One or more atom type properties is invalid in {new_atom_types}.")

        self.atom_types.extend(new_atom_types)


    def _append_topo_types(self, topo_name : str, new_topo_types : List[TopologyType]) -> None:
        """Append new topology types (e.g. types of bond). Skips any that share a name with existing topology."""

        try:
            if len(new_topo_types) == 0:
                logger.warning(f"Just appended a new {topo_name} type of length 0!")
                return
        except:
            raise TypeError(f"new_topo_types must be a list or other iterable, not {type(new_topo_types)}.")
        
        # If the list contains anything that isn't an AtomType object, throw an error.
        if False in [isinstance(t,TopologyType) for t in new_topo_types]:
            raise TypeError(f"One or more topology type properties is invalid in {new_topo_types}.")
        
        param_list = self._get_topo_type_list(topo_name)

        existing_names = [t.name for t in param_list]

        # Don't add new topology if one already exists with the same name.
        new_topo_types = [t for t in new_topo_types if t.name is not None and t.name not in existing_names]

        param_list.extend(new_topo_types)

    def append_bond_types(self, new_bond_types : List[TopologyType]) -> None:
        self._append_topo_types("bond", new_bond_types)
    
    def append_angle_types(self, new_angle_types : List[TopologyType]) -> None:
        self._append_topo_types("angle", new_angle_types)
    
    def append_dihedral_types(self, new_dihedral_types : List[TopologyType]) -> None:
        self._append_topo_types("dihedral", new_dihedral_types)
    
    def append_improper_types(self, new_improper_types : List[TopologyType]) -> None:
        self._append_topo_types("improper", new_improper_types)

    def add_atom_type(self, mass:float, name:str = None, index : Union[str,int] = "append"):
        """Add or insert a new atom type. If index is not specified, appends the new value."""

        # If we already have an AtomType with this name, skip and warn.
        if name is not None and name in (t.name for t in self.atom_types):
            logger.warning(f"Atom Type with name {name} already exists!")
            return
        
        aType = AtomType(mass=mass, name=name)

        if index == "append" or index == len(self.atom_types):
            self.atom_types.append(aType)
            return
        index = int(index)
        if index < len(self.n_atom_types):
            self.atom_types.insert(index,aType)
            # We just changed all atom type indices above index type_id, so we need to shift all references to them.
            self._shift_atom_types(range(index,self.n_atom_types),shift=+1)
        else:
            raise ValueError(f"type_id must be between 0 and the number of existing atom types ({self.n_atom_types}), inclusive. Received {index}.")

    def _add_topo_type(self, topo_name, *params, topo_id = "append", type_name:str=None):
        """
        Add a bond, angle, ... type, with given interaction parameters and optional name.

        If topo_id is unspecified, appends new type.

        If topo_id is an integer, inserts topo type at given index.

        type_name is the name alias of the new bond type.
        """

        # Get list of relevant topology types
        type_list = self._get_topo_type_list(topo_name)

        if type_name is not None and type_name in (t.name for t in type_list):
            logger.warning(f"{topo_name.capitalize()} type with name {type_name} already exists!")
            return

        topoType = TopologyType(*params, name=type_name)

        n_types = len(type_list)

        # Check topo_id has valid value
        if topo_id == "append":
            topo_id = n_types
        elif topo_id not in range(0, n_types+1):
            raise ValueError(f"{topo_name.capitalize()} type index {topo_id} is out of bounds.")
        
        # Offset topology type references in topology list
        self._shift_topo_types(topo_name, range(topo_id,n_types), shift = +1)

        # Insert bond type at relevant index
        type_list.insert(topo_id, topoType)
        

    def add_bond_type(self, *params, bond_id = "append", type_name:str=None):
        self._add_topo_type("bond", *params, topo_id = bond_id, type_name=type_name)
    
    def add_angle_type(self, *params, angle_id = "append", type_name:str=None):
        self._add_topo_type("angle", *params, topo_id = angle_id, type_name=type_name)

    def add_dihedral_type(self, *params, dihedral_id = "append", type_name:str=None):
        self._add_topo_type("dihedral", *params, topo_id = dihedral_id, type_name=type_name)
    
    def add_improper_type(self, *params, improper_id = "append", type_name:str=None):
        self._add_topo_type("improper", *params, topo_id = improper_id, type_name=type_name)

    def _topo_type_exists(self:'World', topo_name:str, type_name:str) -> bool:
        """Does a topology type of name type_name already exist?
        Arguments
        ---------
        topo_name : str
            Topology name (usually bond, angle, dihedral, improper)
        type_name : str
            Topology type identification name
            
        Returns
        -------
        Bool"""

        # If the topology type doesn't have a name, there's no point checking.
        if type_name is None: return False

        type_list = self._get_topo_type_list(topo_name)

        return (type_name in (t.name for t in type_list))
    
    def _validate_topo_atoms(self, *atom_ids):
        """Used to check bond/angle/dihedral atom inputs are valid and order them properly"""
        # Check atom indices are integers
        try:
            atom_ids = [int(a) for a in atom_ids]
        except ValueError as err:
            raise Exception(f"One or more atom indices in {atom_ids} is not an integer.") from err

        # Ensure first atom index is smaller than last
        if atom_ids[0] > atom_ids[-1]:
            atom_ids = atom_ids[::-1]  # Reverses tuple of atoms

        for at in atom_ids:
            if at >= self.n_atoms:
                raise IndexError(f"Atom index {at} is out of bounds.")

        return atom_ids

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
        for topo_list in [self._get_topo_list(topo_name) for topo_name in self._available_topo_types]:
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

    def _get_topo_type_list(self:'World', topo_name) -> List[TopologyType]:
        """Get all stored parameters for a given topology type (e.g. bond, angle)."""
        if not isinstance(topo_name,str):
            raise TypeError(f"Topo name must be a string, not {type(topo_name)}.")

        if topo_name == "bond": return self.bond_types
        elif topo_name == "angle": return self.angle_types
        elif topo_name == "dihedral": return self.dihedral_types
        elif topo_name == "improper": return self.improper_types
        else: raise ValueError(f"{repr(topo_name)} is not a valid value of topo_name.")
 

def lammps_index(i : int):
    """Input a LAMMPS-style (i.e. starting at 1) index and return
        a python-style (0-start) index. Really just takes 1 from an int."""
    
    return i - 1

def _exactly_one_defined(*args):
    """Provide any number of variables,
        returns True if exactly one of them is not None,
        otherwise returns False"""
    return (sum((a is not None for a in args)) == 1)