# This script prints the final expression for the effective average elastic moduli, C,
# according to the Mori-Tanaka model for fiber-reinforced composites.

# The symbols in the equation represent the following physical quantities:
# C:  The effective average elastic moduli tensor of the composite.
# Cm: The fourth-order elasticity tensor of the polymer matrix.
# Cf: The fourth-order elasticity tensor of the fiber.
# Vf: The volume fraction of the fibers.
# Vm: The volume fraction of the matrix (note: Vm = 1 - Vf).
# A:  The Eshelby strain-concentration tensor, which relates the strain
#     in an isolated fiber to the strain in the surrounding matrix.
# I:  The fourth-order identity tensor.
#
# In this notation, multiplication (*) implies a tensor product (e.g., C*A),
# and the exponent **-1 signifies a tensor inversion.

print("C = Cm + Vf * (Cf - Cm) * A * (Vm*I + Vf*A)**-1")