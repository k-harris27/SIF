from SIF import io,system,typing_tools
import argparse

"""
When doing something like crosslinking in MD, we tend to gain extra atom,
    bond, angle, dihedral etc. types over the equilibrated pre-reaction system.
    Rather than messing around rerunning simulations just to have appropriate
    atom types, we can take our equilibrated monomer system and update the atom/
    topo types to match a reference system (which we can use the bond/react template
    creation source data for!).

This script takes (a lot of) arguments to work! Use python topo_extension.py -h for info.
"""

def main():
    
    # Read in arguments from bash.
    args = read_args()

    # Read paths to input, reference and output files.
    d_in = args.data_in
    d_ref = args.data_ref
    d_out = args.data_out
    s_in = args.settings_in
    s_ref = args.settings_ref
    atom_types = read_new_atom_types(args)

    # Read in topology & settings (Forcefield) info of world to have atom/topo types added.
    world = io.read_lammps_data(d_in,s_in)
    
    # Read in topology & settings (Forcefield) info of world used to take atom/topo types from (reference).
    ref = io.read_lammps_data(d_ref,s_ref)
    
    # Atom types have to be added manually currently...
    for t in atom_types:
        world.add_atom_type(type_id=t[0],mass=t[1])
    
    # Update topology types of world to match ref, then checks that there are no actual bonds of the new types.
    # The print lines will break if the topology list breaks - things seem to work fine so can probably
    #   delete them.
    typing_tools.update_types_to_match(world,ref, debug = True)
   
    # Finally, write type-extended data out
    io.write_lammps_data(world, d_out)

def read_args():
    """Interpret the arguments passed to the script when it was run."""
    parser = argparse.ArgumentParser(
		description=("Update the atom/topo types of an input LAMMPS system"
			"to match a reference system.")
		)
    parser.add_argument("-i", "--data-in", action="store", required=True,
			help="Path to input data file.",
			metavar="FILE")
    parser.add_argument("-r", "--data-ref", action="store", required=True,
			help="Path to topology reference data file.",
			metavar="FILE")
    parser.add_argument("-o", "--data-out", action="store", required=False,
			help=("Path to output data file."
				"If not specified, overwrites input file."),
			metavar="FILE")
    parser.add_argument("-I", "--settings-in", action="store", required=True,
			help="Path to input .in.settings file.",
			metavar="FILE")
    parser.add_argument("-R", "--settings-ref", action="store", required=True,
			help="Path to topology reference .in.settings file.",
			metavar="FILE")
    parser.add_argument("-O", "--settings-out", action="store", required=False,
			help=("Path to output .in.settings file."
				"If not specified, overwrites input file."
				"NOTE: Will be identical to ref, so not important."),
			metavar="FILE")
    parser.add_argument("-a", "--atom-type", action="extend", required=False,
			default=None, nargs="+",
			help=("Atom type index (integer) to insert "
				"& corresponding mass (float)."),
			metavar="'TYPE,MASS'")

    return parser.parse_args()


def read_new_atom_types(args):
    """Interpret the atom type list from strings to numbers."""
    list_in = args.atom_type    
    list_out = []

    # If no atom types were assigned, return an empty list.
    if list_in is None: return []
    
    # If atom types were assigned, we need to interpret the strings into a proper list.
    for itm in list_in:
        type_id,mass = itm.split(",")

        # Checking we have sensible values.
        try:
            type_id = int(type_id)
        except ValueError as e:
            raise ValueError(f"Atom type IDs must be integers. "
			     "A non-integer has been received: {type_id}") from e
        try:
            mass = float(mass)
        except ValueError as e:
            raise ValueError(f"Atom mass must be a decimal number. "
			     "A non-decimal has been received: {mass}") from e

        list_out.append((type_id,mass))  

    return list_out

main()
