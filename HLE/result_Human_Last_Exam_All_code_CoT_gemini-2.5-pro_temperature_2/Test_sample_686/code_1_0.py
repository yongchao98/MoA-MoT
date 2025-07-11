def print_susceptibility_formula():
    """
    Prints the derived formula for the magnetic susceptibility chi.
    The formula is expressed in terms of the model parameters and the
    abbreviation N provided in the problem description.
    """

    # Symbolic representation of the variables
    N = "N"
    J = "J"
    beta = "beta"
    c = "c"
    m_cav_sq = "m_cav^2"  # Represents the square of the cavity magnetization

    # Numerator of the main fraction: (c - 1) * J * beta * (1 - m_cav^2)
    numerator = f"({c} - 1) * {J} * {beta} * (1 - {m_cav_sq})"
    
    # Denominator of the main fraction: 1 - (c - 1) * J * beta * (1 - m_cav^2)
    denominator = f"1 - ({c} - 1) * {J} * {beta} * (1 - {m_cav_sq})"

    # Final expression for chi
    chi_expression = f"chi = {N} * (({numerator}) / ({denominator}))"

    # Printing the results
    print("The derived magnetic susceptibility is:")
    print(chi_expression)
    print("")
    print("Each symbol in the equation stands for:")
    print("chi: The magnetic susceptibility")
    print(f"N:   A constant defined as beta * c * (1 - m_0^2) / ({c} - 1)")
    print("c:   The connectivity of the sparse random graph")
    print("J:   The homogeneous coupling constant")
    print("beta: The inverse temperature (1/kT)")
    print("m_0^2: The square of the spontaneous magnetization at a site")
    print("m_cav^2: The square of the spontaneous cavity magnetization, which satisfies the self-consistency equation m_cav = tanh(beta*(c-1)*J*m_cav)")

if __name__ == '__main__':
    print_susceptibility_formula()
