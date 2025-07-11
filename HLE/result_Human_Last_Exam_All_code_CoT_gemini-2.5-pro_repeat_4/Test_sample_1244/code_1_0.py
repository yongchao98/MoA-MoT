import numpy as np

# This script verifies the existence of a vector for question (b).
# It follows the constructive proof outlined above.
# We construct a vector x in Z^14 and check if its squared norm
# is a multiple of 18. This corresponds to the condition
# x_L.x_L is a multiple of 6 for the associated lattice vector x_L.
# The construction ensures the vector is 3-primitive.

# A hypothetical codeword c of weight 6. Any self-dual ternary code
# of length 14 is known to have codewords of weight 6.
c = np.array([1, 1, 1, 1, 1, 1] + [0]*8, dtype=int)

# A vector q chosen to satisfy the necessary congruence conditions.
# For c_0=1, we choose q = 2*e_0.
# This gives q.q = 4 (even) and q.c = 2 (which is 2 mod 3).
q = np.array([2] + [0]*13, dtype=int)

# The vector x in Z^14 is constructed as x = 3q + c.
x = 3 * q + c

# Calculate the squared norm of x.
x_dot_x = np.dot(x, x)

# The condition on the lattice vector norm is x_L.x_L = 0 (mod 6).
# This is equivalent to (x.x)/3 = 0 (mod 6), which means x.x = 0 (mod 18).
# The vector is 3-primitive because c is not the zero vector.
print(f"For question (b), we check our constructed vector x.")
print(f"The equation for our vector is x = 3*q + c")
print(f"q = {q}")
print(f"c = {c}")
print(f"x = {x}")
print(f"The final equation for the squared norm is x.x = {x_dot_x}.")
print(f"The result {x_dot_x} is divisible by 18, so the condition is satisfied.")
