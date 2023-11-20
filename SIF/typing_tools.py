from copy import copy
from .io import *
from ._core import *
from typing import List,Union

def update_types_to_match(target : World, reference : World, debug : bool = False) -> None:
    """Update the atom and topology types of target to be the same as reference.\n
    Apply to a system where all types in target are in reference, but not all types in reference are in target.\n
    Assumes all topology already in target is ordered the same as in the reference world.\n
    WARNING: Atom types can only be updated if all atoms have associated type labels."""

    # Atoms - only if all atom types in both worlds have names
    all_type_lists = [w._get_topo_type_list(topo_kind) for w in (target,reference) for topo_kind in w._available_topo_types]
    all_type_lists += [target.atom_types, reference.atom_types]
    if any(t.name is None for types in all_type_lists for t in types):
        raise ValueError("Some types do not currently have type labels! Type merging can't be confidently carried out without these.")
    names_in_target = [t.name for t in target.atom_types]
    for i,atom_type in enumerate(reference.atom_types):
        if atom_type.name not in names_in_target:
            logger.debug(f"Inserting atom type {atom_type.name} to position {i}.")
            target.add_atom_type(atom_type.mass, atom_type.name, index=i)

    # Topology
    for topo_kind in World._available_topo_types:
        targ_types = target._get_topo_type_list(topo_kind)
        ref_types = reference._get_topo_type_list(topo_kind)

        names_in_target = [t.name for t in target._get_topo_type_list(topo_kind)]
        for i,topo_type in enumerate(ref_types):
            if topo_type.name in names_in_target:
                if topo_type.parameters != targ_types[names_in_target.index(topo_type.name)].parameters:
                    logger.warning(f"{topo_kind} types with the same name but different parameters were found during type updating! Using target world values.")
            else:
                logger.debug(f"Inserting {topo_kind} params {topo_type} at ID {i}.")
                target._add_topo_type(topo_kind, *topo_type.parameters, topo_id = i, type_name=topo_type.name)  # Insert new topology type & update all existing topology to the correct types.

def change_atom_type(world : World, new_type : Union[int,str,AtomType], target_type : Union[int,str] = None, target_indices : Iterable[int] = None,):
    """
    Change the types of atoms matching a range of different requirements. Directly modifies the supplied world.\n
    NOTE: Has no effect on any topology - this is mainly intended for being able to distinguish atoms in post-processing.

    Parameters
    ----------
    world : World
        World object containing atoms etc.
    new_type : int | str | AtomType
        Atom type integer/label to set selected atoms to. Can be 0<=i<n_atom_types, or a new AtomType object, which is appended if given.
    target_type : int | string (optional)
        Select atoms with this original atom type, defined by other type ID (int) or type name (string).
    target_indices : Iterable[int] (optional)
        Select atoms with indices in this iterable only.

    Returns
    -------
    None
    """
    
    if not isinstance(world,World):
        raise TypeError(f"world object must be of type World, received type {type(world)}.")
    
    if not isinstance(new_type,int):
        raise TypeError(f"new_type must be an integer, received {type(new_type)}.")
    
    if not isinstance(target_type,int):
        raise TypeError(f"target_type must be an integer, received {type(target_type)}.")
    

    if isinstance(new_type,AtomType):
        world.add_atom_type(new_type)
    elif new_type >= world.n_atom_types or new_type < 0:
        raise ValueError("new_type index out of bounds! Make sure you're using indices starting from 0.")
    
    # By default, we'll iterate over all atoms
    atoms = world.atoms
    # If we've set specific atom indices to look at, we'll only iterate over those.
    if target_indices is not None:
        atoms = [world.atoms[i] for i in target_indices]

    target_type_given = target_type is not None

    # If an atom type name string was given, find the associated type ID number.
    if isinstance(target_type,str):
        target_type = world.atom_type_id_from_name(target_type)

    for atom in atoms:
        # If a target type is set but doesn't match current atom type,
        #   that is the only case the type is *not* changed. 
        if target_type_given and atom.type_id != target_type:
             continue  # Skip this atom.
        
        # In all other cases, we reassign the atom type id.
        atom.type_id = new_type
            
    return
