import math

# The problem boils down to finding the number of boolean functions of two
# variables that are true for exactly half of the possible inputs.
# This is a combinatorial problem of choosing k items from a set of n.

# n is the total number of rows in the truth table for two variables (2^2).
n = 4
# k is the number of 'True' outcomes required for the probability to be 1/2.
k = 2

# To fulfill the requirement of outputting each number in the final equation,
# we first calculate the factorial components of the combination formula C(n, k).
# The formula is C(n, k) = n! / (k! * (n-k)!)
n_factorial = math.factorial(n)
k_factorial = math.factorial(k)
n_minus_k_factorial = math.factorial(n - k)

# The final result is the combination C(4, 2).
result = math.comb(n, k)

# Print the final equation with all its numerical components.
print(f"The number of fixed points corresponds to the number of ways to choose k=2 'True' states from n=4 total states.")
print(f"The equation is C(n, k) = n! / (k! * (n-k)!).")
print(f"Substituting the numbers:")
print(f"C({n}, {k}) = {n}! / ({k}! * ({n}-{k})!)")
print(f"C({n}, {k}) = {n_factorial} / ({k_factorial} * {n_minus_k_factorial})")
print(f"C({n}, {k}) = {n_factorial} / ({k_factorial * n_minus_k_factorial})")
print(f"C({n}, {k}) = {result}")
print(f"\nThus, there are {result} Quantum Temporal Fixed Points.")
