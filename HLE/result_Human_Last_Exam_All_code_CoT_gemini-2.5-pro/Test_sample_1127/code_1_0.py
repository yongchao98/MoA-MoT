import numpy as np

# We aim to find the minimal polynomial of the connective constant, mu.
# We interpret the connective constant as the spectral radius of the graph's
# adjacency operator. This is derived from the eigenvalues of the
# k-space adjacency matrix A(k).
# The eigenvalues are given by the functions:
# lambda = 2*c +/- sqrt(2 + 2*c), where c = cos(k) ranges in [-1, 1].

# We define the two eigenvalue functions.
def f1(c):
    return 2*c + np.sqrt(2 + 2*c)

def f2(c):
    return 2*c - np.sqrt(2 + 2*c)

# We find the maximum absolute value of these functions over the interval c in [-1, 1].
# From mathematical analysis:
# The range of f1(c) on [-1, 1] is [-2, 4].
# The range of f2(c) on [-1, 1] is [-9/4, 0].
# The complete spectrum is the union of these ranges: [-9/4, 4].
# The spectral radius is the maximum absolute value in the spectrum.
mu_val = max(abs(-9.0/4.0), abs(4.0))

# The calculated value of the connective constant is an integer.
# mu = 4.
# For a rational number q, the minimal polynomial over Q is x - q.
# Therefore, the minimal polynomial is x - 4.
# The equation is x - 4 = 0.

print(f"The calculated connective constant is mu = {int(mu_val)}.")
print("The minimal polynomial P(x) of mu is x - 4.")
print("The equation P(mu) = 0 is expressed with its coefficients as:")

# The polynomial is 1*x - 4 = 0.
coeff_x = 1
coeff_const = -4
print(f"{coeff_x} * {int(mu_val)} + ({coeff_const}) = 0")
