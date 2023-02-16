from argparse import ArgumentError
from typing import List,Set
from .io import *
from ._core import *

def get_connected_atoms(world : World, atom_ids : List[int], extent : int = 3, return_equivalences : bool = False) -> World:
    """
    Returns a new world containing only the atoms connected by (extent) number of bonds
    to the atoms of ids given in atom_ids.

    Arguments
    ---------
    world : World
        The World containing the atoms of interest
    atom_ids : int or List[int]
        The atom id(s) being used to test proximity
    extent : int (default = 3)
        The maximum number of bonds an atom can be from the reference atoms to be included
    return_equivalences : bool
        If true, a dictionary containing the conversions from all atoms to fragment atoms IDs is returned also.
    """

    if not isinstance(world,World):
        raise TypeError(f"world must be a World object, not {type(world)}.")

    if isinstance(atom_ids,int):
        # Test if atom_ids is an int. If it is, put it in a 1-element list
        atom_ids = [atom_ids]
    else:
        try:
            # Maybe it's a list, see if it is and if every value is an int-type
            atom_ids = [int(i) for i in atom_ids]
        except:
            raise TypeError(f"Atom IDs must be either an int or list of ints, not {type(atom_ids)}.")
        
    try:
        extent = int(extent)
    except:
        raise TypeError(f"extent must be an int, not {type(extent)}.")

    keep_atoms = set(atom_ids)
    for atom in atom_ids:
        keep_atoms.update(depth_limited_search(world,depth=3,start_atom=atom))

    return world.get_isolated_atoms(keep_atoms,return_equivalences=return_equivalences)
    

def equate_charges(
        product : World,
        charge : float = float('nan'),
        reactants : World = None,
        dp : int = 5,
        const_atoms : list = []
        ) -> None:
    """Modify the charges of product atoms so that the total charge is equal
    to either charge (float) or the total charge of reactants (World).

    Directly modifies the charges of product, so returns None.
    
    Arguments
    ---------
    product : World
        World containing exclusively fragment atoms needing charge-averaging.
    charge : float
        Target sum total charge (optional).
    reactants : World
        World containing reactive fragment atoms, used to calculate target charge (optional).
    dp : int
        Number of decimal places of accuracy for resulting charges.
    const_atoms : list<int>
        Atom indices of atoms to keep at constant charge, i.e. excluded from offsetting."""
    
    # ^ = xor, error only raised if neither or both are defined.
    if math.isnan(charge) ^ bool(reactants):
        raise ArgumentError("Only one of charge or reactants arguments can be specified.")

    # If we are given a world containing the reactant fragments, calculate their total charge.
    if reactants:
        charge = reactants.total_charge
    
    delta_charge = product.total_charge - charge  # Amount total charge needs to change by

    # Calculate number of atoms we're spreading the charge offset across.
    num_prod_atoms = product.n_atoms - len(const_atoms)
    offset = -delta_charge / num_prod_atoms  # Calculate charge offset per atom.

    # Instantiate cascade rounder to allow rounding without changing the total.
    cascade_rounder = CascadeRound(dp = 5)

    # Offset charges & round to required dp without changing total (cascade rounding).
    for i, atom in enumerate(product.atoms):
        if i in const_atoms: continue  # Skip atoms marked as being kept constant charge.
        atom.charge += offset
        atom.charge = cascade_rounder.round(atom.charge)

    # If the charges aren't equal, raise an error.
    prod_charge = product.total_charge # Debug
    prod_atom_charges = [a.charge for a in product.atoms] # Debug
    assert product.total_charge == charge, "A charge difference still exists between the products and reactants! Something is wrong..."

    del cascade_rounder
    return


def depth_limited_search(world : World, depth : int, start_atom : int, prev_atom : int = -1) -> Set[int]:
    """Recursive function for depth limited search of atoms within a given number of bonds. Don't set prev_atom - only used internally."""

    # Initialisation
    #if prev_atom == -1: depth_limited_search.found = {start_atom}  # Curly braces make a Set instead of List

    found = {start_atom}

    # If we've reached the max distance, don't do anything, just return current atom
    if depth <= 0: return found

    for a in world.find_bonded_atoms(start_atom):
        if a != prev_atom and a not in found:
            found.update(depth_limited_search(world, depth-1,start_atom = a, prev_atom = start_atom))

    return found


class CascadeRound:
    """Object used to keep track of an instance of cascade rounding.\n
    Cascade rounding ensures total remains constant while rounding numbers."""
    def __init__(self, dp : int):
        self.total_f = 0  # Pre-rounding total
        self.total_r = 0  # Rounded total

        # Check dp is an integer variable.
        assert isinstance(dp,int), f"Decimal places must be an int, not {type(dp)}"
        self.dp = dp
    
    def round(self, val : float):
        self.total_f += val
        rounded = round(self.total_f - self.total_r, self.dp)
        self.total_r += rounded
        return rounded
