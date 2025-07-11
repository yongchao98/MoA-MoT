def solve_ising_susceptibility():
    """
    This function prints the derived formula for the magnetic susceptibility chi.
    The formula is expressed in terms of the constant N, the connectivity c,
    and the propagation factor T.
    """
    
    # Define the symbols used in the equation as strings
    N = "N"
    c = "c"
    T = "T"
    
    # Construct the final equation string
    # The numbers 1 are explicitly part of the formula
    numerator = f"{N} * ({c} - 1) * {T}"
    denominator = f"1 - ({c} - 1) * {T}"
    
    # Print the final result
    print(f"chi = ({numerator}) / ({denominator})")
    
    # Print the definitions of the constants
    print("\nwhere:")
    print(f"N = beta * c * (1 - m_0**2) / (c - 1)")
    print(f"T = tanh(beta*J) * (1 - m_cav**2) / (1 - m_msg**2)")

solve_ising_susceptibility()