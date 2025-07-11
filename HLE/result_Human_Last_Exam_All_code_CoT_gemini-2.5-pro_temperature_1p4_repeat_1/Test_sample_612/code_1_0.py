# This script prints the expression for the effective average elastic moduli (C)
# for a fiber-reinforced composite based on the Mori-Tanaka model.
# The formula is presented symbolically.

# In the context of tensor algebra:
# - Products are tensor contractions.
# - The inverse of a tensor T is denoted as inv(T).

print("Based on the Mori-Tanaka model, the expression for the effective average elastic moduli tensor C is:")
print("")
# The following line prints each component of the final equation.
print("C = Cm + Vf * (Cf - Cm) * A * inv(Vm * I + Vf * A)")
print("")
print("Where the terms are:")
print("C:  Effective average elastic moduli tensor of the composite")
print("Cm: Elasticity tensor of the polymer matrix")
print("Cf: Elasticity tensor of the fiber")
print("Vm: Volume fraction of the matrix")
print("Vf: Volume fraction of the fiber")
print("A:  Eshelby strain-concentration tensor")
print("I:  Fourth-order identity tensor")