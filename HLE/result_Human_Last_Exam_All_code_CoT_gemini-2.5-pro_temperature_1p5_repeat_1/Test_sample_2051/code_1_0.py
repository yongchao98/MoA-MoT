def calculate_polynomial_degree(g, n_plus, n_minus):
    """
    Calculates the degree of the volume polynomial Z_{g, n_+, n_-}.

    Args:
        g (int): The genus of the surface.
        n_plus (int): The number of positively oriented boundaries.
        n_minus (int): The number of negatively oriented boundaries.
    """
    
    # The total number of boundaries n
    n = n_plus + n_minus
    
    # The degree d' of the polynomial in terms of L_i^2 is 3g - 3 + n.
    # This formula gives the dimension of the moduli space of curves of genus g with n marked points.
    degree_in_L_squared = 3 * g - 3 + n
    
    # The degree of the polynomial in terms of L_i is twice the degree in terms of L_i^2.
    final_degree = 2 * degree_in_L_squared
    
    # --- Output the step-by-step calculation ---
    print(f"For the case g = {g}, n_+ = {n_plus}, n_- = {n_minus}:")
    print(f"1. The total number of boundaries is n = n_+ + n_- = {n_plus} + {n_minus} = {n}.")
    
    print("\n2. The volume is a polynomial in the squared lengths L_i^2.")
    print(f"   Its degree, d', is given by the formula: d' = 3*g - 3 + n.")
    print(f"   Substituting the values: d' = 3*{g} - 3 + {n} = {degree_in_L_squared}.")

    print("\n3. To find the degree in terms of the lengths L_i, we use d = 2 * d'.")
    print(f"   d = 2 * {degree_in_L_squared} = {final_degree}.")
    
    print(f"\nConclusion: The degree of the polynomial Z_{{{g},{n_plus},{n_minus}}}(L_+ | L_-) is {final_degree}.")

# --- Parameters from the question (b) ---
g_val = 0
n_plus_val = 3
n_minus_val = 1

# --- Run the calculation and print the result ---
calculate_polynomial_degree(g_val, n_plus_val, n_minus_val)
