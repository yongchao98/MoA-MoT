# Define the symbols used in the Mori-Tanaka model expression
C = "C"        # Effective average elastic moduli of the composite
Cf = "Cf"      # Fourth-order elasticity tensor of the fiber
Cm = "Cm"      # Fourth-order elasticity tensor of the polymer matrix
Vf = "Vf"      # Volume fraction of the fibers
Vm = "Vm"      # Volume fraction of the matrix
A = "A"        # Eshelby strain-concentration tensor
I = "I"        # Fourth-order identity tensor

# Construct the expression for C using the Mori-Tanaka model.
# The ":" symbol denotes a tensor double-dot product (contraction).
# The "inv(...)" function denotes tensor inversion.
# The "*" symbol denotes scalar multiplication.
equation = f"{C} = {Cm} + {Vf} * ({Cf} - {Cm}) : {A} : inv({Vm} * {I} + {Vf} * {A})"

# Print the final expression
print("The Mori–Tanaka model expression for the effective average elastic moduli (C) is:")
print(equation)

# Print the definition of each term for clarity
print("\nWhere the symbols in the equation represent:")
print(f"  {C}: The effective average elastic moduli of the composite")
print(f"  {Cm}: The fourth-order elasticity tensor of the polymer matrix")
print(f"  {Cf}: The fourth-order elasticity tensor of the fiber")
print(f"  {Vm}: The volume fraction of the matrix (Vm = 1 - Vf)")
print(f"  {Vf}: The volume fraction of the fibers")
print(f"  {A}: The Eshelby strain-concentration tensor")
print(f"  {I}: The fourth-order identity tensor")
