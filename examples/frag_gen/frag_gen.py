from utils import io,system,template_tools
import math

"""
In this script, we use SIF to load in a LAMMPS data file containing
    a number of polymer building blocks to build fix bond/react
    template files from. A number of commands need to be used to
    properly create templates, so we go through each and explain
    its purpose.
"""

# Read in the world containing all minimum n-mers from data file.
all_mols = io.read_lammps_data("poly_atoms.data")

# We have identified the desired "reacting" atoms of each molecule ahead of time
#   so we make a note of their (0-starting) atom index in the data file
#   (and corresponding index within the individual molecule) here.
# Epoxy reactive atoms (epoxide O = 11 & primary C = 10)
# Linker reactive atoms (Nitrogen = 49(0), Amine H = 59(10))
# DGEBA-MXDA dimer reactive atoms (Primary C = 82(11), Nitrogen = 83(12))
# DGEBA-DGEBA-MXDA trimer "reactive" atoms
#   (Primary C = 153(11), Nitrogen = 154(12), Primary C = 164(22))


# The example simulations contain dihedral interactions, so we need to
#   include all atoms within 3 bonds of the reacting atoms.
#   get_connected_atoms() takes a source world (all_mols), list of reactive atom
#   IDs in the source world, and "extent" (the maximum number of bonds from the
#   reactive atoms an atom can be to be included in the template).

# Monomers - Useful for calculations later but not used directly
dgeba_frag = template_tools.get_connected_atoms(all_mols,[10, 11], extent=3)
mxda_frag = template_tools.get_connected_atoms(all_mols,[49, 59], extent=3)

# Reactant Fragments - Include both reactant fragments for intermolecular reactions.
#   The equiv_atoms dictionaries aren't used in this example, but is useful if you
#   want to be able to easily get from fragment atom IDs to the source world IDs.
dgeba_mxda_frag, dgeba_mxda_equiv_atoms = template_tools.get_connected_atoms(
        all_mols,[10, 11, 49, 59], extent=3, return_equivalences=True)
dgeba_dm_frag, dgeba_dm_equiv_atoms = template_tools.get_connected_atoms(
        all_mols,[10,11,82, 83], extent=3, return_equivalences=True)

# Products - Make sure you get the exact same parts of the molecules as the reactants.
dm_frag, dm_equiv_atoms = template_tools.get_connected_atoms(
        all_mols,[82, 83], extent=3, return_equivalences=True)
ddm_frag, ddm_equiv_atoms = template_tools.get_connected_atoms(
        all_mols,[153,154,164], extent=3, return_equivalences=True)


### Charge equalisation between reactants and products ###
#
# Since there is an Aromatic C overlapping between MXDA templates on either
#   side, we have to keep that atom's charge consistent with unreacted MXDA.
#   Also, since that atom is indistinguishable from another within the range
#   of the template file, we keep the charges of both constant.
#   To do this, we assign the original monomer charges to the relevant atoms,
#   flag them to be left alone in the charge equilibration step, and add all
#   atoms other than them to a "fragment" within the .template file. We can 
#   then tell fix bond/react to ONLY update the charges of the atoms in the
#   fragment we define. Deep breaths.
###


# Template-indexed atom IDs to keep charge contant
# Currently needs prior knowledge of the fragment atom IDs,
#   but can probably be done without using equivalence dictionaries.
dgeba_mxda_const_q_atoms = [13,14]  # MXDA Aromatic Carbons
dgeba_dm_const_q_atoms = [18,19]    # MXDA Aromatic Carbons
dm_const_q_atoms = [8,9]            # MXDA Aromatic Carbons
ddm_const_q_atoms = [8,9]           # MXDA Aromatic Carbons

# Get the original monomer charges of the atoms (in order wrt above lists).
const_q_charges = [dgeba_mxda_frag.atoms[i].charge for i in dgeba_mxda_const_q_atoms]

# Assign original charges to all constant charge atoms in post-reaction fragments.
# (Mostly to avoid confusion when looking at template atom charges later).
for q,dgeba_dm_i,dm_i,ddm_i in zip(
        const_q_charges,dgeba_dm_const_q_atoms,dm_const_q_atoms,ddm_const_q_atoms):

    dgeba_dm_frag.atoms[dgeba_dm_i].charge = q
    dm_frag.atoms[dm_i].charge = q
    ddm_frag.atoms[ddm_i].charge = q


# Monomer fragment charge totals
dgeba_frag_charge = dgeba_frag.total_charge
mxda_frag_charge = mxda_frag.total_charge

# Calculate fragment target charges
dm_q_target = dgeba_frag_charge + mxda_frag_charge
ddm_q_target = 2* dgeba_frag_charge + mxda_frag_charge

# Fragment charge equalisation.
# Reactant charges aren't used in fix bond/react, so don't matter.
# We tell the code to ignore the atoms whose charges should be kept constant throughout
#   the simulation by giving their template atom IDs to const_atoms.
template_tools.equate_charges(dm_frag, charge = dm_q_target, const_atoms = dm_const_q_atoms)
template_tools.equate_charges(ddm_frag, charge = ddm_q_target, const_atoms = ddm_const_q_atoms)

# Create sub-fragments in reactant templates to define which atoms' charges **are**
# updated in fix bond/react (see custom_charges in fix bond/react docs).
dgeba_mxda_frag.add_fragment(
        "update_q",
        [i for i in range(dgeba_mxda_frag.count_atoms()) if i not in dgeba_mxda_const_q_atoms])
dgeba_dm_frag.add_fragment(
        "update_q",
        [i for i in range(dgeba_dm_frag.count_atoms()) if i not in dgeba_dm_const_q_atoms])


# In the OPLS-AA Forcefield, hydrogen atoms have no Lennard-Jones interactions.
#   This leads to problems with the weird geometries made when crosslink topology is
#   applied, and simulations have a tendency to explode due to hydrogen atoms collapsing
#   onto negative atoms they were previously bonded to. To prevent this, we prevent the H
#   atoms from "seeing" the previously bonded negative atoms by adding a ghost bond
#   with 0 force, since atoms within 2 bonds of each other ignore long-range interactions.
# Again, needs prior knowledge of template atom IDs currently, but could probably be done
#   using the equivalency lists.

dm_n_bond_types = len(dm_frag.bond_params)
dm_frag.add_bond_type(0.,0.)
dm_frag.add_bond(3,5,type_id = dm_n_bond_types)

dgeba_dm_n_bond_types = len(dgeba_dm_frag.bond_params)
dgeba_dm_frag.add_bond_type(0.,0.)
dgeba_dm_frag.add_bond(13,15,type_id=dgeba_dm_n_bond_types)
ddm_n_bond_types = len(ddm_frag.bond_params)

# Only the new alcohol Ox gets a ghost bond, old ghost bond is removed.
# We don't need the stabilisation it provides any more, so may as well.
# Be careful when assigning equivalent atoms in the .map file later though!
# The ghost bond needs to be on the newly bonded atoms, not the previously
#   bonded atoms.
# (see REACTER introduction.pptx notes.)
ddm_frag.add_bond_type(0.,0.)
ddm_frag.add_bond(3,5,type_id = ddm_n_bond_types)

# Output in LAMMPS data file format
io.write_lammps_data(dgeba_mxda_frag,
        "dgeba-mxda_frag.data",
        "DGEBA & MXDA fragment created by SIF.")
io.write_lammps_data(dm_frag,
        "dm_frag.data",
        "DGEBA-MXDA Dimer fragment created by SIF.")
io.write_lammps_data(dgeba_dm_frag,
        "dgeba-dm_frag.data",
        "DGEBA-MXDA Dimer & DGEBA fragment created by SIF.")
io.write_lammps_data(ddm_frag,
        "ddm_frag.data",
        "DGEBA-DGEBA-MXDA Trimer fragment created by SIF.")

# Output in LAMMPS molecule (bond/react template) format
io.write_react_template(dgeba_mxda_frag,
        "test_data/q_eq_frags/dgeba-mxda_frag.template",
        "DGEBA & MXDA fragment created by SIF.")
io.write_react_template(dm_frag,
        "test_data/q_eq_frags/dm_frag.template",
        "DGEBA-MXDA Dimer fragment created by SIF.")
io.write_react_template(dgeba_dm_frag,
        "test_data/q_eq_frags/dgeba-dm_frag.template",
        "DGEBA-MXDA Dimer & DGEBA fragment created by SIF.")
io.write_react_template(ddm_frag,
        "test_data/q_eq_frags/ddm_frag.template",
        "DGEBA-DGEBA-MXDA Trimer fragment created by SIF.")
