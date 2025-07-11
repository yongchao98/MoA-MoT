def solve():
    """
    This function prints the derivation and the final formula for the magnetic susceptibility chi.
    """
    # Define symbols as strings for printing the formula
    N = "N"
    c = "c"
    K = "\uD835\uDCDA" # Script K for the propagation factor
    chi = "\u03C7"    # chi
    
    # Expression for chi in terms of N, c, and K
    numerator = f"({c} - 1) * {K}"
    denominator = f"1 - ({c} - 1) * {K}"
    
    # The final expression for chi
    final_chi_expr = f"{N} * ({numerator}) / ({denominator})"

    # Definitions of the components
    m0 = "m_0"
    beta = "\u03B2"
    T = "tanh({beta}*J)"
    m_cav = "m_cav"
    
    # Print the equation for chi
    print(f"The magnetic susceptibility {chi} is given by the expression:")
    print(f"{chi} = {final_chi_expr}")
    print()

    # Print definitions of the terms
    print("Where:")
    print(f"N = {beta} * c * (1 - {m0}^2) / ({c} - 1)")
    print(f"{K} = {T} * (1 - {m_cav}^2) / (1 - ({T})^2 * {m_cav}^2)")
    print(f"T = tanh({beta}*J) is the hyperbolic tangent of the coupling constant J times beta.")
    print(f"m_cav is the cavity magnetization, which is related to the total magnetization m_0 by the self-consistency equation:")
    print(f"{m0} * (1 + {T} * {m_cav}^2) = {m_cav} * (1 + {T})")

solve()
