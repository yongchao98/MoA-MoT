import math

# The problem is to find the value of foo(7).
# Based on mathematical analysis, foo(n) counts the number of labeled groups of order n.
# For n=7, this is given by the formula: 7! / |Aut(Z_7)|.

# The order of the group.
n = 7

# For a prime number n, the size of the automorphism group of the cyclic group Z_n
# is phi(n) = n - 1.
aut_group_size = n - 1

# Calculate n factorial.
n_factorial = math.factorial(n)

# The result is n! divided by the size of the automorphism group.
result = n_factorial // aut_group_size

# Print the numbers used in the final equation.
# The equation is: 7! / |Aut(Z_7)| = result
# which is equivalent to: 5040 / 6 = 840
print(f"{n_factorial} / {aut_group_size} = {result}")
