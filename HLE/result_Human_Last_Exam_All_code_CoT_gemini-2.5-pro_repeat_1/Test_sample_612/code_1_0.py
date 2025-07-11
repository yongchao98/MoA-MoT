# Define the quantities involved in the Mori-Tanaka model as string variables
# to construct the final expression.
I = "I"   # Fourth-order identity tensor
Cf = "Cf" # Fourth-order elasticity tensor of the fiber
Cm = "Cm" # Fourth-order elasticity tensor of the matrix
Vf = "Vf" # Volume fraction of the fibers
Vm = "Vm" # Volume fraction of the matrix (Note: Vm + Vf = 1)
A = "A"   # Eshelby strain-concentration tensor

# The Mori-Tanaka model provides an expression for the effective average
# elastic moduli, C. The derivation relates the average stress <sigma> to the
# average strain <epsilon> through the properties of the constituent phases.
#
# The final expression is C = (Vm*Cm + Vf*Cf*A) : [Vm*I + Vf*A]^(-1)
# where ':' denotes a tensor contraction (like a product) and '^(-1)' denotes the inverse.

# The following code constructs and prints this formula.
# We will use '*' to represent the tensor contraction for clarity in the output.
print("The Mori-Tanaka expression for the effective elastic moduli C is:")

# Constructing the first part of the expression: (Vm * Cm + Vf * Cf * A)
term1 = f"({Vm} * {Cm} + {Vf} * {Cf} * {A})"

# Constructing the second part (the inverse term): (Vm * I + Vf * A)**(-1)
term2_inv = f"({Vm} * {I} + {Vf} * {A})**(-1)"

# Printing the final equation C = term1 * term2_inv
# This fulfills the requirement to output each symbol in the final equation.
print(f"C = {term1} * {term2_inv}")