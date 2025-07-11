def solve_force_equation():
    """
    This function prints the symbolic equation for the x-directed total force
    on the specified conducting material.
    The derivation is based on calculating the self-force on block 2.
    """
    
    # Define the symbols used in the equation as strings
    force = "F_x"
    equal_sign = "="
    term1 = "-a*D"
    term2 = "mu_0 / 2"
    term3 = "(I_0**2 / D**2)"
    term4 = "(sigma_2 / (sigma_1 + sigma_2))**2"
    
    # Assemble the final equation string
    # F_x = -aD * (mu_0/2) * (I_0^2/D^2) * (sigma_2/(sigma_1+sigma_2))^2
    final_equation = f"{force} {equal_sign} {term1} * ({term2}) * {term3} * {term4}"

    # Print the equation, reconstructing it piece by piece for clarity.
    print(f"The x-directed total force on the conducting material is:")
    print(f"F_x = -aD * (mu_0 / 2) * (I_0^2 / D^2) * (sigma_2 / (sigma_1 + sigma_2))^2")

solve_force_equation()