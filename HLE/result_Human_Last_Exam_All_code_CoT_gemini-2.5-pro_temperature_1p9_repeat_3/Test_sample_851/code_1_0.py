import collections

# Define the properties of all 32 crystallographic point groups.
# 'chiral': Lacks any roto-inversion axis (m, -1, -3, -4, -6).
# 'polar': Has a unique polar direction.
# 'g_tensor_nonzero': Can exhibit gyrotropy (non-zero gyration tensor). This requires being 
# non-centrosymmetric, but some non-centrosymmetric classes like -6 and -6m2 are excluded 
# because their specific symmetries still force the tensor to be zero.
CrystalClass = collections.namedtuple('CrystalClass', ['name', 'chiral', 'polar', 'g_tensor_nonzero'])

point_groups = [
    # Triclinic
    CrystalClass(name="1", chiral=True, polar=True, g_tensor_nonzero=True),
    CrystalClass(name="-1", chiral=False, polar=False, g_tensor_nonzero=False),
    # Monoclinic
    CrystalClass(name="2", chiral=True, polar=True, g_tensor_nonzero=True),
    CrystalClass(name="m", chiral=False, polar=True, g_tensor_nonzero=True),
    CrystalClass(name="2/m", chiral=False, polar=False, g_tensor_nonzero=False),
    # Orthorhombic
    CrystalClass(name="222", chiral=True, polar=False, g_tensor_nonzero=True),
    CrystalClass(name="mm2", chiral=False, polar=True, g_tensor_nonzero=True),
    CrystalClass(name="mmm", chiral=False, polar=False, g_tensor_nonzero=False),
    # Tetragonal
    CrystalClass(name="4", chiral=True, polar=True, g_tensor_nonzero=True),
    CrystalClass(name="-4", chiral=False, polar=False, g_tensor_nonzero=True),
    CrystalClass(name="4/m", chiral=False, polar=False, g_tensor_nonzero=False),
    CrystalClass(name="422", chiral=True, polar=False, g_tensor_nonzero=True),
    CrystalClass(name="4mm", chiral=False, polar=True, g_tensor_nonzero=True),
    CrystalClass(name="-42m", chiral=False, polar=False, g_tensor_nonzero=True),
    CrystalClass(name="4/mmm", chiral=False, polar=False, g_tensor_nonzero=False),
    # Trigonal
    CrystalClass(name="3", chiral=True, polar=True, g_tensor_nonzero=True),
    CrystalClass(name="-3", chiral=False, polar=False, g_tensor_nonzero=False),
    CrystalClass(name="32", chiral=True, polar=False, g_tensor_nonzero=True),
    CrystalClass(name="3m", chiral=False, polar=True, g_tensor_nonzero=True),
    CrystalClass(name="-3m", chiral=False, polar=False, g_tensor_nonzero=False),
    # Hexagonal
    CrystalClass(name="6", chiral=True, polar=True, g_tensor_nonzero=True),
    CrystalClass(name="-6", chiral=False, polar=False, g_tensor_nonzero=False),
    CrystalClass(name="6/m", chiral=False, polar=False, g_tensor_nonzero=False),
    CrystalClass(name="622", chiral=True, polar=False, g_tensor_nonzero=True),
    CrystalClass(name="6mm", chiral=False, polar=True, g_tensor_nonzero=True),
    CrystalClass(name="-6m2", chiral=False, polar=False, g_tensor_nonzero=False),
    CrystalClass(name="6/mmm", chiral=False, polar=False, g_tensor_nonzero=False),
    # Cubic
    CrystalClass(name="23", chiral=True, polar=False, g_tensor_nonzero=True),
    CrystalClass(name="m-3", chiral=False, polar=False, g_tensor_nonzero=False),
    CrystalClass(name="432", chiral=True, polar=False, g_tensor_nonzero=True),
    CrystalClass(name="-43m", chiral=False, polar=False, g_tensor_nonzero=True),
    CrystalClass(name="m-3m", chiral=False, polar=False, g_tensor_nonzero=False),
]

# Find classes that are achiral, non-polar, and have a non-zero gyration tensor.
valid_classes = []
for pg in point_groups:
    is_achiral = not pg.chiral
    is_non_polar = not pg.polar
    allows_gyrotropy = pg.g_tensor_nonzero
    
    if is_achiral and is_non_polar and allows_gyrotropy:
        valid_classes.append(pg.name)

print("The crystal classes that are achiral, non-polar, and have the correct symmetry for gyrotropic effects are:")
for cls in sorted(valid_classes):
    print(cls)
