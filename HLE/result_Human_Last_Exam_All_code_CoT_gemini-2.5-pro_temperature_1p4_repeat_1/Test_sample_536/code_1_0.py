import math

# Let A be the inner product of h_p and b_p, and B be the inner product of h_p and z_p.
A = 0.9375
B = 0.9

# Let L be the inner product of b_p and z_p that we want to find.
# The coplanarity of the three vectors h, b, z implies that the determinant of their Gram matrix is zero.
# det([[1, A, B], [A, 1, L], [B, L, 1]]) = 0
# This expands to the quadratic equation: L^2 - 2*A*B*L + (A^2 + B^2 - 1) = 0.

# The solutions for L are given by the quadratic formula:
# L = A*B +/- sqrt((A*B)^2 - (A^2 + B^2 - 1))
# L = A*B +/- sqrt(A^2*B^2 - A^2 - B^2 + 1)
# L = A*B +/- sqrt((1 - A^2)*(1 - B^2))

# This corresponds to the geometric relation L = cos(arccos(A) +/- arccos(B)).
# A reasonable physical assumption is that the vector h lies on the great circle arc
# connecting b and z on the unit sphere. This implies that the angle between b and z
# is the sum of the angles between (b, h) and (h, z).
# This selects the solution L = cos(arccos(A) + arccos(B)) = AB - sqrt((1-A^2)(1-B^2)).

A_squared = A**2
B_squared = B**2
term1 = 2 * A * B
term2 = A_squared + B_squared - 1

# Calculate the chosen root for L
L = A * B - math.sqrt((1 - A_squared) * (1 - B_squared))

# We print the equation with the calculated numerical coefficients
# The equation is L^2 - (term1)*L + (term2) = 0
print(f"The relationship between the inner products leads to the quadratic equation for L = <b_p, z_p>:")
print(f"L^2 - (2 * {A} * {B}) * L + ({A**2} + {B**2} - 1) = 0")
print(f"L^2 - {term1} * L + {term2} = 0")
print("\nSolving for L, we get:")
print(f"L = {L}")
