import math

# The number of variables in the Boolean function
n = 4

# For a(n), k is the floor of n/2
k = math.floor(n / 2)

# Calculate the components of the binomial coefficient formula C(n, k) = n! / (k! * (n-k)!)
n_factorial = math.factorial(n)
k_factorial = math.factorial(k)
n_minus_k_val = n - k
n_minus_k_factorial = math.factorial(n_minus_k_val)

# Calculate the final result
result = n_factorial // (k_factorial * n_minus_k_factorial)

# Print the step-by-step calculation
print(f"The maximal number of prime implicants, a(n), for a Boolean function of n variables is C(n, floor(n/2)).")
print(f"For n = {n}, we need to calculate a({n}) = C({n}, {k}).\n")
print("The calculation is as follows:")
print(f"a({n}) = C({n}, {k}) = {n}! / ({k}! * ({n}-{k})!)")
print(f"a({n}) = {n_factorial} / ({k_factorial} * {n_minus_k_factorial})")
print(f"a({n}) = {n_factorial} / ({k_factorial * n_minus_k_factorial})")
print(f"a({n}) = {result}")
