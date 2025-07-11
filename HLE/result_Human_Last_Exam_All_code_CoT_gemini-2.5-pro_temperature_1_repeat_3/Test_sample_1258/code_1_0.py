import sys

def find_demagnetizing_factor_expression():
    """
    This function constructs and prints the analytical expression for the
    fluxmetric demagnetizing factor of a cylinder based on the given parameters.
    """

    # Define the symbols used in the expression based on the problem description
    N_f = "N_f"  # Fluxmetric demagnetizing factor
    g = "g"      # Length-to-diameter ratio
    k = "k"      # Modulus for elliptic integrals
    F_k = "F(k)" # Complete elliptic integral of the first kind
    E_k = "E(k)" # Complete elliptic integral of the second kind

    # Construct the terms of the expression as strings
    # The numbers in the formula are 8, 2, 1, and 4.
    # Term 1: (8 / (g^2 * k^2)) * [E(k) - (1 - k^2) * F(k)]
    num_8 = "8"
    g_sq = f"{g}^2"
    k_sq = f"{k}^2"
    term1_coeff = f"({num_8} / ({g_sq} * {k_sq}))"
    one_minus_k_sq = f"(1 - {k_sq})"
    elliptic_term = f"[{E_k} - {one_minus_k_sq} * {F_k}]"
    term1 = f"{term1_coeff} * {elliptic_term}"

    # Term 2: (4 / g^2)
    num_4 = "4"
    term2 = f"({num_4} / {g_sq})"

    # Combine the terms to form the final expression
    final_expression = f"{N_f} = {term1} - {term2}"

    # Print the results
    print("The analytical expression for the fluxmetric demagnetizing factor (N_f) for a magnetic cylinder is:")
    print(final_expression)
    print("\nwhere:")
    print(f"  {g} is the length-to-diameter ratio (L/D).")
    print(f"  {k} is the modulus for the elliptic integrals, defined as {k}^2 = 1 / (1 + {g}^2 / 4).")
    print(f"  {F_k} is the complete elliptic integral of the first kind with modulus {k}.")
    print(f"  {E_k} is the complete elliptic integral of the second kind with modulus {k}.")

if __name__ == '__main__':
    find_demagnetizing_factor_expression()