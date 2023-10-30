import SIF

sys = SIF.World(-10,-10,-10,10,10,10)

# Ideally have actual atom types to add.
[sys.add_atom_type(SIF.AtomType(1.00, f"H{i}")) for i in range(3)]

[print(t.name) for t in sys.atom_types]

sys.add_atom(1,-1.4,[-1.,-3.,4.])
sys.add_atom(2,-1.4,[-2,-3.,4.])
sys.add_atom(0,1.4,[0,0,0])
sys.add_atom(1,1.0,[0,0,0])


sys.add_bond_type(0.0, 0.0, type_name="H0-H0")

[print(b.name) for b in sys.bond_types]

sys.add_bond(0,2, "H0-H0")
sys.add_bond(0,1, "H0-H0")
sys.add_bond(2,3, "H0-H0")

sys.add_angle_type(0.0, 0.0)
sys.add_angle_type(0.0, 0.0)

[print(a.name) for a in sys.angle_types]

sys.add_angle(3,1,2, type_id=0)
sys.add_angle(3,0,2, type_id=1)
sys.add_angle(2,0,3, type_id=1)

sys.delete_atom(1)

atom = sys.atoms[0].copy()

print(sys.bonds)

print(sys.angles)

del(sys)
