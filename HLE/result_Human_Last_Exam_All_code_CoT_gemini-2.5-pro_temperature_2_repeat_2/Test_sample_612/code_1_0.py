# The Mori-Tanaka model is used to determine the effective elastic moduli of a composite material.
# The following symbols represent the given quantities:
# C:  The effective average elastic moduli tensor of the composite
# I:  The fourth-order identity tensor
# Cf: The fourth-order elasticity tensor of the fiber
# Cm: The fourth-order elasticity tensor of the polymer matrix
# Vf: The volume fraction of the fibers
# Vm: The volume fraction of the matrix (Vm = 1 - Vf)
# A:  The Eshelby strain-concentration tensor

# The model derives the effective tensor C by relating the average stress and strain in the composite.
# The resulting expression for C is printed below.
# In the expression, '*' denotes tensor contraction (multiplication) and 'inv()' represents the inverse of a tensor.

print("C = Cm + Vf * (Cf - Cm) * A * inv(Vm * I + Vf * A)")