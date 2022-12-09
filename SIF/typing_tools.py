from copy import copy
from .io import *
from .system import *

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
