import math

def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code
    on a sphere with n smooth and m rough holes.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    if n < 0 or m < 0 or not isinstance(n, int) or not isinstance(m, int):
        print("Error: n and m must be non-negative integers.")
        return

    # The number of logical qubits k is given by n + m - 2
    # This formula is based on the simplified model described.
    k = n + m - 2
    
    print(f"For n = {n} smooth holes and m = {m} rough holes:")
    
    # We need to handle the case where k might be negative for small n, m
    if k < 0:
        print(f"The formula k = n + m - 2 = {n} + {m} - 2 = {k} results in a negative number of qubits.")
        print("This indicates the simplified model is not applicable for these small values.")
        print("The actual GSD for such cases is typically 1 or 2, corresponding to 0 or 1 qubits.")
    else:
        # GSD = 2^k
        gsd = int(math.pow(2, k))
        print(f"The number of logical qubits is k = n + m - 2 = {n} + {m} - 2 = {k}.")
        print(f"The ground space degeneracy is 2^k = 2^{k} = {gsd}.")

# Example usage with n=4 smooth holes and m=3 rough holes.
n_smooth = 4
m_rough = 3
calculate_toric_code_gsd(n_smooth, m_rough)
