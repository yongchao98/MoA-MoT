def get_potential_distribution_expression():
    """
    This function prints the derived expression for the Electrical Double-Layer (EDL)
    potential distribution psi(y) based on the problem's conditions.
    """
    
    # Define the variables symbolically for printing the equation
    psi_y = "psi(y)"
    z1 = "z_1"
    beta = "beta"
    k = "k"
    H = "H"
    y = "y"
    
    # Construct the parts of the final equation as strings
    zeta_1_expr = f"{z1} * (1 + {beta} * {k})"
    sinh_term = f"sinh({k} * ({H} - {y}))"
    denominator = f"sinh({k} * {H})"
    
    # Assemble the final expression
    final_expression = f"{psi_y} = ( {zeta_1_expr} * {sinh_term} ) / ( {denominator} )"
    
    # Print the explanation and the final answer
    print("The expression for the Electrical Double-Layer potential distribution psi(y) is:")
    print(final_expression)
    print("\nWhere:")
    print(f"  {psi_y}: The electrical potential at a distance y from the bottom surface.")
    print(f"  {z1}: A parameter related to the zeta potential of the bottom surface.")
    print(f"  {beta}: The slip length.")
    print(f"  {k}: The Debye-Huckel parameter.")
    print(f"  {H}: The height of the microchannel.")
    print(f"  {y}: The vertical position within the channel (0 <= y <= H).")
    print(f"  sinh: The hyperbolic sine function.")
    print(f"  The number '1' is included in the term for the slip-dependent zeta potential.")

# Execute the function to display the result
get_potential_distribution_expression()

# The final expression is psi(y) = ( z_1 * (1 + beta * k) * sinh(k * (H - y)) ) / ( sinh(k * H) )
final_answer_string = "psi(y) = (z_1*(1 + beta*k) * sinh(k*(H - y))) / sinh(k*H)"
# For the purpose of providing a single answer string.
# print(f"<<<{final_answer_string}>>>")