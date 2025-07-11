def calculate_toric_code_degeneracy(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code on a torus
    with n smooth holes and m rough holes.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n < 0 or m < 0:
        print("Error: The number of holes (n and m) must be non-negative integers.")
        return

    # The number of logical qubits k is n + m.
    # The ground space degeneracy is 2^k.
    exponent = n + m
    degeneracy = 2**exponent

    # Print the equation with the specific numbers, as requested.
    print(f"The ground space degeneracy is 2^({n} + {m}) = {degeneracy}")

# Example usage with n=2 smooth holes and m=3 rough holes.
n_smooth_holes = 2
m_rough_holes = 3
calculate_toric_code_degeneracy(n_smooth_holes, m_rough_holes)