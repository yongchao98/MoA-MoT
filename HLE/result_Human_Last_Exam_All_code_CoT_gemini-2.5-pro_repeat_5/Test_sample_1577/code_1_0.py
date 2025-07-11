def calculate_toric_code_degeneracy(n, m):
    """
    Calculates and explains the ground state degeneracy (GSD) of the toric code
    on a surface with n smooth holes and m rough holes.

    This calculation assumes the surface has genus g=0 and that n>=1 and m>=1.

    Args:
        n (int): The number of smooth holes (boundaries).
        m (int): The number of rough holes (boundaries).
    """
    if n <= 0 or m <= 0:
        print("This formula GSD = 2^(n+m-2) is derived for the case where n>=1 and m>=1.")
        return

    # The number of logical qubits, k, is the sum of qubits from
    # smooth boundaries (n-1) and rough boundaries (m-1).
    k = n + m - 2

    # The ground state degeneracy is 2 to the power of k.
    gsd = 2**k

    # Print the explanation and the final equation with the numbers substituted.
    print(f"Given n = {n} smooth holes and m = {m} rough holes:")
    print("The number of logical qubits (k) is calculated as: k = (n - 1) + (m - 1)")
    print(f"k = ({n} - 1) + ({m} - 1) = {n + m - 2}")
    print(f"The ground state degeneracy (GSD) is 2^k.")
    print(f"GSD = 2^({n} + {m} - 2) = 2^{k} = {gsd}")

# --- User-configurable values ---
# Set the number of smooth and rough holes here
num_smooth_holes = 4
num_rough_holes = 3
# --- End of configuration ---

# Execute the calculation
calculate_toric_code_degeneracy(num_smooth_holes, num_rough_holes)
