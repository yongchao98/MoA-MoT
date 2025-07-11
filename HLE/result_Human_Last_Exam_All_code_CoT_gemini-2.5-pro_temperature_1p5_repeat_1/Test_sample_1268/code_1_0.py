import math

# This script states the relationship between the upper bound for an ideal's
# norm in any ideal class (the Minkowski bound M_K) and the covolume V
# for a real quadratic field K = Q(sqrt(N)).

# The upper bound is proportional to the covolume V.
# The constant of proportionality is derived from the Minkowski bound formula
# for a number field of degree n=2 with r2=0 pairs of complex embeddings.
# The formula for the constant is: C = (n! / n^n)

# For a real quadratic field, the degree is n=2.
n = 2
constant_factor = math.factorial(n) / (n**n)

# The relationship is: k_k,inf <= M_K = C * V
print("The upper bound for the maximum norm (k_k,inf) in relation to the covolume (V) is given by the Minkowski bound M_K.")
print("The relationship is of the form: k_k,inf <= C * V")
print("\nFor a real quadratic field, the degree is n=2. The constant C is calculated as:")
print(f"C = n! / n^n = {n}! / {n}^{n}")

print("\nTherefore, the final equation showing the upper bound is:")
# The instruction is to output each number in the final equation.
# The number is the constant factor 0.5.
print(f"k_k,inf <= {constant_factor} * V")