def solve():
    """
    Calculates the largest value k such that for any valid arrangement of k diamonds,
    at least one diamond can be moved to an adjacent cell without violating the rules.
    This is equivalent to M-1, where M is the size of the smallest "frozen" configuration.
    
    A minimal "frozen" configuration on an n x n grid (n=2m, even) can be constructed
    by taking a base grid of diamonds on all (2i, 2j) cells and adding patches on
    the borders to prevent moves.
    
    The size of the base grid (S_ee) is m * m.
    The number of patches needed for the right border is m/2.
    The number of patches needed for the bottom border is m/2.
    
    Total size M = m^2 + m/2 + m/2 = m^2 + m.
    The required value k is M - 1.
    """
    n = 2024
    m = n // 2
    
    # Size of the base configuration S_ee
    size_S_ee = m * m
    
    # Number of patch diamonds needed for one edge (right or bottom)
    num_patches_per_edge = m // 2
    
    # Size of the minimal frozen configuration M
    M = size_S_ee + num_patches_per_edge + num_patches_per_edge
    
    # The result k is M - 1
    k = M - 1
    
    # The final equation is m*m + m - 1
    print(f"The grid size is {n}x{n}.")
    print(f"We can denote m = n/2 = {m}.")
    print("The size of the smallest frozen arrangement is M = m^2 + m.")
    print(f"M = {m}^2 + {m} = {m*m} + {m} = {M}")
    print("The largest value k such that any arrangement is not frozen is M - 1.")
    print(f"k = {M} - 1 = {k}")

solve()