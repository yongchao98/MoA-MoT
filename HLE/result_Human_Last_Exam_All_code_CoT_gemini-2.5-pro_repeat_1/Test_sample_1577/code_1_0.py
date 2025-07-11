import math

def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code
    on a planar surface with n smooth holes and m rough holes.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    # The number of logical qubits k is n + m for a planar surface
    # with n smooth and m rough holes (which implies n+m+1 total boundaries).
    # The GSD is 2^k.
    
    k = n + m
    gsd = 2**k

    print(f"For a toric code with:")
    print(f"  n = {n} (smooth holes)")
    print(f"  m = {m} (rough holes)")
    print("-" * 25)
    print("The ground space degeneracy (GSD) is calculated as 2^(n + m).")
    print(f"GSD = 2^({n} + {m})")
    print(f"GSD = 2^{k}")
    print(f"GSD = {gsd}")


# Example values for n and m.
# You can change these to explore other cases.
n_smooth_holes = 2
m_rough_holes = 1

calculate_toric_code_gsd(n_smooth_holes, m_rough_holes)
