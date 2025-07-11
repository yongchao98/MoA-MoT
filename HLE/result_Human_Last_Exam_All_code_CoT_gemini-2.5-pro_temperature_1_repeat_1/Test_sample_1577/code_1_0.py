# Number of smooth holes
n = 3

# Number of rough holes
m = 4

# The ground space degeneracy of the toric code on a sphere with n smooth and m rough holes
# is given by 2^k, where k is the number of logical qubits.
# Assuming n >= 1 and m >= 1, the number of logical qubits k is (n-1) + (m-1) = n+m-2.
# This corresponds to answer choice B.

if n >= 1 and m >= 1:
    # Calculate the exponent k
    k = n + m - 2
    
    # Calculate the ground space degeneracy
    degeneracy = 2**k
    
    # Print the equation step-by-step
    # Remember in the final code you still need to output each number in the final equation!
    print(f"For n = {n} smooth holes and m = {m} rough holes:")
    print(f"The ground space degeneracy is 2^(n + m - 2) = 2^({n} + {m} - 2) = 2^{k} = {degeneracy}")
else:
    # The more general formula is k = max(0, n-1) + max(0, m-1)
    k_general = max(0, n - 1) + max(0, m - 1)
    degeneracy_general = 2**k_general
    print(f"The formula 2^(n+m-2) assumes n>=1 and m>=1.")
    print(f"Using the general formula for n={n}, m={m}:")
    print(f"k = max(0, {n}-1) + max(0, {m}-1) = {k_general}")
    print(f"Degeneracy = 2^{k_general} = {degeneracy_general}")
