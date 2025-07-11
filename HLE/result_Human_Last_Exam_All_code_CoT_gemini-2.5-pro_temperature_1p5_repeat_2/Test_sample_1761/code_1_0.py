import math

# The dimension 'n' of the complex projective space P^n.
# You can change this value to any integer greater than 0.
n = 5

# The complex dimension we want to compute is h^0(P^n, Omega^1(2)).
# From the Euler sequence, we derive the formula:
# dim = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))
# where h^0(P^n, O(k)) = C(n+k, k), with C being the binomial coefficient "n choose k".

# First term: dimension of H^0(P^n, O(1)^(n+1))
# This is equal to (n+1) * h^0(P^n, O(1))
# h^0(P^n, O(1)) = C(n+1, 1)
h0_O1 = math.comb(n + 1, 1)
dim_term1 = (n + 1) * h0_O1

# Second term: dimension of H^0(P^n, O(2))
# h^0(P^n, O(2)) = C(n+2, 2)
dim_term2 = math.comb(n + 2, 2)

# The final dimension is the difference between the two terms.
result = dim_term1 - dim_term2

# Print the final equation with all the computed numbers.
print(f"For the projective space P^n with n = {n}, the calculation is as follows:")
print(f"Dimension of H^0(P^{n}, O(1)^(n+1)) = (n+1) * C(n+1, 1) = {n+1} * {h0_O1} = {dim_term1}")
print(f"Dimension of H^0(P^{n}, O(2)) = C(n+2, 2) = C({n+2}, 2) = {dim_term2}")
print("\nThe final dimension is the difference between these two values.")
print(f"Final equation: {dim_term1} - {dim_term2} = {result}")