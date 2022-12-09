from ..utils import system

sys = system.World(-10,-10,-10,10,10,10)

sys.add_atom(1,-1.4,[-1.,-3.,4.])
sys.add_atom(2,-1.4,[-2,-3.,4.])
sys.add_atom(0,1.4,[0,0,0])
sys.add_atom(1,1.0,[0,0,0])

sys.add_bond_type(0.0, 0.0)

sys.add_bond(0,2)
sys.add_bond(0,1)
sys.add_bond(2,3)

sys.add_angle_type(0.0, 0.0)
sys.add_angle_type(0.0, 0.0)

sys.add_angle(3,1,2)
sys.add_angle(3,0,2, type_id=1)
sys.add_angle(2,0,3, type_id=1)

sys.delete_atom(1)

atom = sys.atoms[0].copy()

print(sys.bonds)

print(sys.angles)

del(sys)