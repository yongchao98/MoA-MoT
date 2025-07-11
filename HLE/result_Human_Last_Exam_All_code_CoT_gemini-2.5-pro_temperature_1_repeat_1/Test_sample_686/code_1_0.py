def print_susceptibility_formula():
    """
    Prints the derived formula for the magnetic susceptibility chi.
    """
    # Define the components of the formula as string variables
    beta = "beta"
    c = "c"
    m0_sq = "m_0^2"
    K = "K"
    J = "J"
    m_cav_sq = "m_cav^2"
    
    # Construct the numerator and denominator of the final expression
    # Numerator: beta * c * (1 - m_0^2) * K
    numerator_str = f"{beta} * {c} * (1 - {m0_sq}) * {K}"
    # Denominator: 1 - (c - 1) * K
    denominator_str = f"1 - ({c} - 1) * {K}"
    
    # Print the final equation and the definitions of the terms
    print("The final expression for the magnetic susceptibility is:")
    print(f"chi = ({numerator_str}) / ({denominator_str})")
    print("\nWhere the terms are defined as:")
    
    # Print each term in the final equation
    print(f"Numerator = {numerator_str}")
    print(f"Denominator = {denominator_str}")
    
    print(f"\nAnd the variables represent:")
    print(f"  {beta}: The inverse temperature (1/kT)")
    print(f"  {c}: The connectivity of the random graph")
    print(f"  m_0: The magnetization at a reference site 0")
    print(f"  {K}: The perturbation propagation factor, defined as:")
    k_numerator = f"tanh({beta}*{J}) * (1 - {m_cav_sq})"
    k_denominator = f"1 - (tanh({beta}*{J}))^2 * {m_cav_sq}"
    print(f"     K = ({k_numerator}) / ({k_denominator})")
    print(f"     J: The homogeneous coupling constant")
    print(f"     m_cav: The cavity magnetization")

# Execute the function to print the result
print_susceptibility_formula()