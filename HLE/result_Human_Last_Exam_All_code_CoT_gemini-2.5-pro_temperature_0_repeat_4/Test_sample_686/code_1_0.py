def solve_ising_susceptibility():
    """
    This function prints the final derived equation for the magnetic susceptibility chi.
    The equation is expressed in terms of N, the connectivity c, and the propagation factor P.
    """
    # Define the symbolic variables in the equation
    chi = "χ"
    N = "N"
    c = "c"
    P = "P"
    
    # The numbers present in the final equation
    one = 1
    
    # Construct the numerator and denominator of the fraction
    # Numerator: N * (c - 1) * P
    # Denominator: 1 - (c - 1) * P
    
    # Print the final equation
    # The format is chi = Numerator / Denominator
    print(f"{chi} = ({N} * ({c} - {one}) * {P}) / ({one} - ({c} - {one}) * {P})")

solve_ising_susceptibility()