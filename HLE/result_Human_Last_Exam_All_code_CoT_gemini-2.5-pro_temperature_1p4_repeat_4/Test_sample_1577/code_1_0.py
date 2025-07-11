import math

def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code on a torus
    with n smooth holes and m rough holes.

    Args:
        n (int): The number of smooth holes (must be non-negative).
        m (int): The number of rough holes (must be non-negative).
    """
    if n < 0 or m < 0:
        print("Error: Number of holes cannot be negative.")
        return

    # The problem implies there are holes, so n + m > 0.
    # On a torus (g=1), the number of logical qubits k is n + m.
    # The ground space degeneracy is 2^k.
    # Let's check the n=0, m=0 case for completeness (a closed torus).
    # In that case, k=2g=2. The formula k=n+m gives 0, so it's not applicable.
    # The problem is about surfaces with holes, so we assume n+m > 0.
    if n + m == 0:
        k = 2 # Special case for a closed torus with no holes
        gsd = 2**k
        print("For a closed torus with n=0 and m=0:")
        print(f"The number of logical qubits is k = 2.")
        print(f"The ground space degeneracy is 2^2 = {gsd}")
    else:
        k = n + m
        # Use try-except for potentially large numbers, though standard floats handle this
        try:
            gsd = 2**k
            # We want to print the full equation as requested
            print(f"For a surface with n={n} smooth holes and m={m} rough holes, the ground space degeneracy is:")
            print(f"2^({n}+{m}) = {gsd}")
        except OverflowError:
            print(f"The result 2^{k} is too large to represent as a standard float.")


# Example usage:
# Let's use n=2 smooth holes and m=1 rough hole as an example.
n_smooth_holes = 2
m_rough_holes = 1
calculate_toric_code_gsd(n_smooth_holes, m_rough_holes)
