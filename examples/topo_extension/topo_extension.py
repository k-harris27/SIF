import SIF
import SIF.io
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

    # Read in topology & settings (Forcefield) info of world to have atom/topo types added.
    world = SIF.io.read_lammps_data(d_in,s_in)
    
    # Read in topology & settings (Forcefield) info of world used to take atom/topo types from (reference).
    ref = SIF.io.read_lammps_data(d_ref,s_ref)
    
    # Update topology types of world to match ref, then checks that there are no actual bonds of the new types.
    # The print lines will break if the topology list breaks - things seem to work fine so can probably
    #   delete them.
    SIF.typing_tools.update_types_to_match(world,ref, debug = True)
   
    # Finally, write type-extended data out
    SIF.io.write_lammps_data(world, d_out)

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
    parser.add_argument("-o", "--data-out", action="store", required=True,
			help=("Path to output data file."),
			metavar="FILE")
    parser.add_argument("-I", "--settings-in", action="store", required=True,
			help="Path to input .in.settings file.",
			metavar="FILE")
    parser.add_argument("-R", "--settings-ref", action="store", required=True,
			help="Path to topology reference .in.settings file.",
			metavar="FILE")
    """
    parser.add_argument("-p", "--no-pair-coeffs", action="store_true",
			help=("If this flag is set, atom pair coefficients "
				"will not be added to the output files."))
    """

    return parser.parse_args()

main()
