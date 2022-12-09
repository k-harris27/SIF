from utils import io,system,template_tools
import math

# World containing all minimum n-mers
all_mols = io.read_lammps_data("test_data/poly_atoms.data")

# Epoxy reactive atoms (epoxide O = 11 & primary C = 10)
# Linker reactive atoms (Nitrogen = 49(0), Amine H = 59(10))
# DGEBA-MXDA dimer reactive atoms (Primary C = 82(11), Nitrogen = 83(12))
# DGEBA-DGEBA-MXDA trimer "reactive" atoms (Primary C = 153(11), Nitrogen = 154(12), Primary C = 164(22))

# Monomers
dgeba_frag = template_tools.get_connected_atoms(all_mols,[10, 11], extent=3)
mxda_frag = template_tools.get_connected_atoms(all_mols,[49, 59], extent=3)

# Reactant Fragments
dgeba_mxda_frag, dgeba_mxda_equiv_atoms = template_tools.get_connected_atoms(all_mols,[10, 11, 49, 59], extent=3, return_equivalences=True)
dgeba_dm_frag, dgeba_dm_equiv_atoms = template_tools.get_connected_atoms(all_mols,[10,11,82, 83], extent=3, return_equivalences=True)

# Products
dm_frag, dm_equiv_atoms = template_tools.get_connected_atoms(all_mols,[82, 83], extent=3, return_equivalences=True)
ddm_frag, ddm_equiv_atoms = template_tools.get_connected_atoms(all_mols,[153,154,164], extent=3, return_equivalences=True)

# Fragment charge assignment sub-fragments (look at custom_charges paragraph in fix bond/react docs)
dgeba_mxda_const_q_atoms = [13,14]  # MXDA Aromatic Carbons
dgeba_dm_const_q_atoms = [18,19]    # MXDA Aromatic Carbons
dm_const_q_atoms = [8,9]            # MXDA Aromatic Carbons
ddm_const_q_atoms = [8,9]           # MXDA Aromatic Carbons
const_q_charges = [dgeba_mxda_frag.atoms[i].charge for i in dgeba_mxda_const_q_atoms]  # In order w/ above

# Assign original charges to all constant charge atoms in post-reaction fragments.
# (Mostly so I don't get confused looking at template atom charges later).
for q,dgeba_dm_i,dm_i,ddm_i in zip(const_q_charges,dgeba_dm_const_q_atoms,dm_const_q_atoms,ddm_const_q_atoms):
    dgeba_dm_frag.atoms[dgeba_dm_i].charge = q
    dm_frag.atoms[dm_i].charge = q
    ddm_frag.atoms[ddm_i].charge = q


### Charge equalisation between reactants and products
# NOTE: Since there is an Aromatic C overlapping between MXDA templates on either
#       side, we have to keep that atom's charge consistent with unreacted MXDA.
###

# Monomer fragment charge totals
# Take away the charges of any atoms to be ignored in charge assignment so we can equilibrate the rest of the charges.
dgeba_frag_charge = dgeba_frag.total_charge
mxda_frag_charge = mxda_frag.total_charge

# Calculate fragment target charges
dm_q_target = dgeba_frag_charge + mxda_frag_charge
ddm_q_target = 2* dgeba_frag_charge + mxda_frag_charge


# Fragment charge equalisation. Reactant charges aren't used, so don't matter.
template_tools.equate_charges(dm_frag, charge = dm_q_target, const_atoms = dm_const_q_atoms)
template_tools.equate_charges(ddm_frag, charge = ddm_q_target, const_atoms = ddm_const_q_atoms)

# Create sub-fragments in reactant fragments to define which atoms' charges **are** updated (see custom_charges in fix bond/react docs).
dgeba_mxda_frag.add_fragment("update_q", [i for i in range(dgeba_mxda_frag.count_atoms()) if i not in dgeba_mxda_const_q_atoms])
dgeba_dm_frag.add_fragment("update_q", [i for i in range(dgeba_dm_frag.count_atoms()) if i not in dgeba_dm_const_q_atoms])

# Add ghost bonds
dm_n_bond_types = len(dm_frag.bond_params)
dm_frag.add_bond_type(0.,0.)
dm_frag.add_bond(3,5,type_id = dm_n_bond_types)

dgeba_dm_n_bond_types = len(dgeba_dm_frag.bond_params)
dgeba_dm_frag.add_bond_type(0.,0.)
dgeba_dm_frag.add_bond(13,15,type_id=dgeba_dm_n_bond_types)
ddm_n_bond_types = len(ddm_frag.bond_params)

# Only the new alcohol Ox gets a ghost bond, old ghost bond is removed - see REACTER introduction.pptx notes.
ddm_frag.add_bond_type(0.,0.)
ddm_frag.add_bond(3,5,type_id = ddm_n_bond_types)

# Output
io.write_lammps_data(dgeba_mxda_frag,"test_data/q_eq_frags/dgeba-mxda_frag.data","DGEBA & MXDA fragment created by SIF.")
io.write_lammps_data(dm_frag,"test_data/q_eq_frags/dm_frag.data", "DGEBA-MXDA Dimer fragment created by SIF.")
io.write_lammps_data(dgeba_dm_frag,"test_data/q_eq_frags/dgeba-dm_frag.data","DGEBA-MXDA Dimer & DGEBA fragment created by SIF.")
io.write_lammps_data(ddm_frag,"test_data/q_eq_frags/ddm_frag.data","DGEBA-DGEBA-MXDA Trimer fragment created by SIF.")

io.write_react_template(dgeba_mxda_frag,"test_data/q_eq_frags/dgeba-mxda_frag.template","DGEBA & MXDA fragment created by SIF.")
io.write_react_template(dm_frag,"test_data/q_eq_frags/dm_frag.template", "DGEBA-MXDA Dimer fragment created by SIF.")
io.write_react_template(dgeba_dm_frag,"test_data/q_eq_frags/dgeba-dm_frag.template","DGEBA-MXDA Dimer & DGEBA fragment created by SIF.")
io.write_react_template(ddm_frag,"test_data/q_eq_frags/ddm_frag.template","DGEBA-DGEBA-MXDA Trimer fragment created by SIF.")