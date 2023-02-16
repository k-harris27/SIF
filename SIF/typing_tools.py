from copy import copy
from .io import *
from ._core import *
from typing import List

def update_types_to_match(target : World, reference : World, debug : bool = False) -> None:
    """Update the ""atom"" and topo types of target to be the same as reference.\n
    Apply to a system where all types in target are in reference, but not all types in reference are in target.\n
    Assumes all topology already in target is ordered the same as in the reference world.\n
    WARNING: Cannot currently differentiate between different atom types of the same mass, so atoms are not handled. Not sure how to???"""

    for topo_name in ("bond","angle","dihedral","improper"):
        targ_params = target._get_topo_param_list(topo_name)
        ref_params = reference._get_topo_param_list(topo_name)

        if len(targ_params) > len(ref_params):
            raise ValueError(f"Reference world ({len(ref_params)} {topo_name}) should have more types than target world ({len(targ_params)} {topo_name}).")

        current_target_param = 0
        targ_params_old = copy(targ_params)  # We can edit targ_params on the go while keeping a copy of the original
        for i,params in enumerate(ref_params):
            if params == targ_params_old[current_target_param]:
                current_target_param += 1
            else:
                if debug: print(f"DEBUG: Inserting {topo_name} params {params} at ID {i}.")
                target._add_topo_type(topo_name, *params, topo_id = i)  # Insert new topology type & update all existing topology to the correct types.

def change_atom_type(world : World, new_type : int, target_type : int = None, target_indices : List[int] = None, new_mass : float = None):
    """
    Change the types of atoms matching a range of different requirements. Directly modifies the supplied world.\n
    NOTE: Has no effect on any topology - this is mainly intended for being able to distinguish atoms in post-processing.

    Parameters
    ----------
    world : World
        World object containing atoms etc.
    new_type : int
        Atom type integer to set selected atoms to. Can be 0<=i<=n_atom_types, with a new atom type being created if new_type=n_atom_types.
    target_type : int (optional)
        Select atoms with this original atom type.
    target_indices : List[int] (optional)
        Select atoms with indices in this list only.
    new_mass : float (optional)
        Mass of newly created atom type. Ignored if an atom type is not created.

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
    
    if new_type > world.n_atom_types:
        raise ValueError("new_type cannot be larger than the number of atom types! Make sure you're using indices starting from 0.")
    elif new_type == world.n_atom_types:
        world.atom_type_params.append(new_mass)
    
    # By default, we'll iterate over all atoms
    atoms = world.atoms
    # If we've set specific atom indices to look at, we'll only iterate over those.
    if target_indices is not None:
        atoms = [world.atoms[i] for i in target_indices]

    for atom in atoms:
        # If a target type is set but doesn't match current atom type,
        #   that is the only case the type is *not* changed. 
        if target_type is not None and atom.type_id != target_type:
             continue  # Skip this atom.
        
        # In all other cases, we reassign the atom type id.
        atom.type_id = new_type
            
    return
