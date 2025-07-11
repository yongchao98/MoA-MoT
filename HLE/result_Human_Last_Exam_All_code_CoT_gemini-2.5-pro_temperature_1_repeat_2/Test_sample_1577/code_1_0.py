def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code
    with n smooth holes and m rough holes on a planar surface.

    Args:
        n (int): The number of smooth holes (must be >= 1).
        m (int): The number of rough holes (must be >= 1).
    """
    if n < 1 or m < 1:
        print("This formula assumes n >= 1 and m >= 1.")
        # Handle cases for n=0 or m=0 if necessary, though the problem implies n,m > 0.
        # The degeneracy is 2^(k-1) for k holes of the same type.
        # This derivation of 2^(n+m-2) is for the mixed-boundary case.
        return

    # The number of logical qubits is (n-1) + (m-1)
    exponent = n + m - 2
    
    # The GSD is 2 to the power of the number of logical qubits
    gsd = 2**exponent

    print(f"The ground space degeneracy for n={n} smooth holes and m={m} rough holes is:")
    print(f"GSD = 2^({n} + {m} - 2) = 2^{exponent} = {gsd}")


# Example usage with n=4 smooth holes and m=3 rough holes.
n_smooth_holes = 4
m_rough_holes = 3
calculate_toric_code_gsd(n_smooth_holes, m_rough_holes)
