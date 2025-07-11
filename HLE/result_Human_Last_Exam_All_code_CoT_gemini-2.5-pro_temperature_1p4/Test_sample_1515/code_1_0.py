import sympy

def print_nsvz_beta_function():
    """
    This function defines and prints the symbolic NSVZ beta function formula.
    It demonstrates the structure for a general N=1 supersymmetric gauge theory
    and then for a specific case of SU(Nc) with Nf flavors in the fundamental rep.
    """

    # Define the symbols used in the general formula
    alpha, T_adj, T_R, gamma, pi, N_f = sympy.symbols('alpha T_adj T_R gamma pi N_f')

    # General Formula: Assume all Nf matter fields are in the same representation R
    # with the same anomalous dimension gamma.
    numerator = 3 * T_adj - N_f * T_R * (1 - gamma)
    denominator = 1 - T_adj * alpha / pi
    beta_alpha = - (alpha**2 / (2 * pi)) * (numerator / denominator)

    print("The symbolic NSVZ Beta Function, beta(alpha), is:")
    print("-------------------------------------------------")
    sympy.pprint(beta_alpha)
    print("\nWhere:")
    print("  alpha = g^2 / (4*pi) is the fine-structure constant.")
    print("  T_adj is the Dynkin index for the adjoint representation (gluons).")
    print("  T_R is the Dynkin index for the matter representation.")
    print("  N_f is the number of matter flavors.")
    print("  gamma is the anomalous dimension of the matter fields.")
    print("  pi is the mathematical constant pi.")

    # Specific case: SU(Nc) gauge theory with Nf flavors in the fundamental representation
    # For SU(Nc): T_adj = Nc, T_R for fundamental = 1/2
    N_c = sympy.symbols('N_c')
    
    # Substitute the values for SU(Nc)
    # The numbers in this final equation are 3, 1/2 (as a fraction), 1, 2, and 1
    beta_su_nc = beta_alpha.subs({T_adj: N_c, T_R: sympy.Rational(1, 2)})

    print("\nFor a specific case of SU(Nc) gauge theory with Nf fundamental flavors:")
    print("T(adj) is replaced with Nc, and T(R) is replaced with 1/2.")
    print("----------------------------------------------------------------------")
    # Print the equation with its terms
    print("beta(alpha) = - (alpha**2 / (2 * pi)) * ( (3 * N_c - N_f * (1/2) * (1 - gamma)) / (1 - N_c * alpha / pi) )")


if __name__ == '__main__':
    print_nsvz_beta_function()
