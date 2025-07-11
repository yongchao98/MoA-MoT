import numpy as np

# Set up a generic problem with m=3, n=2
m, n = 3, 2
# Use a fixed seed for reproducibility
np.random.seed(0)

A = np.random.rand(m, n)
x = np.random.rand(n)
b = np.random.rand(m)

# The simplest candidate for the perturbation E is given by E = r x^T / ||x||^2
# where r = b - Ax. This E makes x an exact solution to (A+E)x = b,
# which is a special case of being a least-squares solution.
# The rank of this E is 1 if r is not the zero vector.

# Calculate the residual r = b - Ax
r = b - np.dot(A, x)

# Check if the residual is the zero vector
is_r_zero = np.allclose(r, np.zeros(m))

print(f"A = \n{A}")
print(f"x = {x}")
print(f"b = {b}")
print(f"Residual r = b - Ax = {r}")

if is_r_zero:
    # This case is highly unlikely for random inputs.
    # If r is zero, E can be the zero matrix, so rank is 0.
    print("The residual is zero. E can be the zero matrix, rank(E) = 0.")
    greatest_rank = 0
else:
    # If r is non-zero, the candidate matrix E = r x^T / ||x||^2 has rank 1.
    # The more general solution for E also turns out to be rank at most 1.
    # Therefore, the greatest possible rank is 1.
    print("The residual is non-zero. A rank-1 solution for E exists.")
    greatest_rank = 1

print(f"\nThe greatest possible rank of E is {greatest_rank}.")