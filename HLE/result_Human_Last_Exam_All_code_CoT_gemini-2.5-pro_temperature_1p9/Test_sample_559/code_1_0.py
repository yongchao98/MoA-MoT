import sympy as sp

# The separatrix for the given system of differential equations is a curve
# in the (d, u) phase plane. Through analysis, this curve can be found
# to satisfy a polynomial equation.
#
# The equation has the form:
# A * u**3 + B * u * d + C * d = 0
#
# The following code defines and prints the coefficients A, B, and C.

# Coefficients of the separatrix equation
A = 3
B = 2
C = 1

# Print the equation in a readable format
u, d = sp.symbols('u d')
equation = A * u**3 + B * u * d + C * d
print("The equation for the separatrix is:")
print(f"{equation} = 0")

# As requested, output each number in the final equation
print("\nThe coefficients of the equation are:")
print(f"Coefficient of u**3 (A): {A}")
print(f"Coefficient of u*d  (B): {B}")
print(f"Coefficient of d    (C): {C}")

# Final Answer Check
# The equation can be written as d*(B*u + C) = -A*u**3
# d = -A*u**3 / (B*u + C)
# Let's substitute the coefficients to get the explicit form for d(u).
# d = -3*u**3 / (2*u + 1)
# We can verify that as u -> 1, d -> -3/(2+1) = -1.
# Also, as u -> 0, d -> 0.
# This matches the saddle point (-1, 1) and the equilibrium (0, 0)
# which the separatrix connects.