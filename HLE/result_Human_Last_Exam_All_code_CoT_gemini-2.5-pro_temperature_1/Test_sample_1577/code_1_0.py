def calculate_toric_code_degeneracy(n, m):
    """
    Calculates and prints the formula for the ground space degeneracy
    of the toric code with n smooth holes and m rough holes.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    if n < 0 or m < 0:
        print("Number of holes cannot be negative.")
        return
    
    # The number of logical qubits k is n + m - 2.
    # The ground space degeneracy is 2^k.
    # This formula is generally valid for n + m >= 2.
    
    print(f"The ground space degeneracy (GSD) of the toric code on a sphere")
    print(f"with n={n} smooth holes and m={m} rough holes is given by the formula: GSD = 2^(n + m - 2).")
    print("\nSubstituting the given values, the equation is:")
    
    # The prompt requires printing each number in the final equation.
    print(f"2^({n}+{m}-2)")

# Example usage with n=5 smooth holes and m=3 rough holes.
n_smooth_holes = 5
m_rough_holes = 3
calculate_toric_code_degeneracy(n_smooth_holes, m_rough_holes)