def solve_susceptibility_formula():
    """
    This function prints the derived formula for the magnetic susceptibility chi.
    """

    # Define the symbols used in the final formula as strings
    var_N = "N"
    var_c = "c"
    var_K = "K"
    var_chi = "χ"
    
    # Numbers that appear in the formula
    num_1_in_c = "1"
    num_1_in_denom = "1"

    # Construct the numerator and denominator strings
    # Numerator: N * (c-1) * K
    numerator = f"{var_N} * ({var_c} - {num_1_in_c}) * {var_K}"
    
    # Denominator: 1 - (c-1) * K
    denominator = f"{num_1_in_denom} - ({var_c} - {num_1_in_c}) * {var_K}"
    
    # Print the final equation
    print("The final expression for the magnetic susceptibility χ is:")
    print(f"{var_chi} = ({numerator}) / ({denominator})")
    
    # Print the definitions of the symbols
    print("\nWhere the terms are defined as:")
    print(f"N = β * c * (1 - m_0^2) / (c - 1)  (as given in the hint)")
    print(f"K = tanh(βJ) * (1 - m_cav^2) / (1 - m_cav^2 * tanh(βJ)^2)")
    print("β is the inverse temperature, J is the coupling constant, c is the connectivity.")
    print("m_0 is the site magnetization and m_cav is the cavity magnetization.")

# Execute the function to print the result
solve_susceptibility_formula()