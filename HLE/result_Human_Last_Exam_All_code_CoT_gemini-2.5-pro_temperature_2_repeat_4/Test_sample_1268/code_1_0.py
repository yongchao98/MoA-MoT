import math

# We choose a valid squarefree natural number N. For a definite quaternion
# algebra over Q, N must be a product of an odd number of primes.
# Let's choose N = 5 (a single prime).
N = 5

# 1. Calculate the covolume (V) of the lattice.
V = N / 2.0

# 2. Calculate the upper bound for the maximum norm (k_k,inf), which is
# its exact value, the covering radius R. The formula is R = sqrt(N)/2.
# Let's call this quantity `k_k_inf`.
k_k_inf = math.sqrt(N) / 2.0

# The relationship is k_k,inf = sqrt(V/2).
# We will now print the final equation with the numbers substituted in.
# Let's present the relationship in its squared form: k_k,inf^2 = V/2

print("For a squarefree natural number N, we analyze the associated lattice.")
print(f"Let's use the example N = {N}.")
print(f"The covolume is V = N / 2 = {N} / 2 = {V}.")
print("")
print("The maximum norm (k_k,inf) is related to the covolume (V) by the equation:")
print("k_k,inf^2 = V / 2")
print("")
print("Substituting the values for N = 5:")
# We output each number in the final equation.
print(f"k_k_inf^2 = {V} / 2 = {V/2}")
print(f"k_k_inf = sqrt({V/2}) = {k_k_inf}")