import math

# The dimension of the complex projective space. n is not specified in the
# problem, so we choose n=3 as a non-trivial example.
n = 3

print(f"Calculating the dimension for the complex projective space P^n where n = {n}.")

# First term in the subtraction: dim H^0(P^n, O(1)^(n+1)) = (n+1)^2
term1 = (n + 1)**2

# Second term in the subtraction: dim H^0(P^n, O(2)) = C(n+2, 2)
term2 = math.comb(n + 2, 2)

# The result is the difference
result = term1 - term2

# Print out the final equation with the computed values
print(f"The dimension is given by the equation: h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))")
print(f"For n={n}, the equation with calculated values is:")
print(f"{term1} - {term2} = {result}")

print(f"\nThe complex dimension of the space of global sections is {result}.")
