import math

# As derived from the algebraic factorization, the four roots of the
# polynomial are sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).

# Create a list with the four roots
roots = [
    math.sqrt(14),
    2 * math.sqrt(6),
    math.sqrt(34),
    2 * math.sqrt(11)
]

# Sort the roots in increasing order
roots.sort()

# The final equation can be expressed in factored form as (X - r1)(X - r2)(X - r3)(X - r4) = 0,
# where r1, r2, r3, and r4 are the roots. The following are these four roots.
print("The four roots of the equation in increasing order are:")
print(roots[0])
print(roots[1])
print(roots[2])
print(roots[3])