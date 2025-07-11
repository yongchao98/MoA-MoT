def print_sliding_angle_equation():
    """
    This function prints the derived equation for the angle at which the
    rod begins to slide, based on the principles of static equilibrium.
    """
    # Define the variables in the final equation as strings for printing.
    var_theta = "theta"
    var_mu = "mu"
    
    # The final equation derived from the force balance is tan(theta) = mu.
    # The prompt mentions outputting numbers in the equation, but this specific
    # result contains only variables. We will present it in its standard, simplified form.
    final_equation = f"tan({var_theta}) = {var_mu}"
    
    print("The expression for the angle theta at which the rod begins to slide is given by the equation:")
    print(final_equation)
    
    print("\nThis can also be written as:")
    solution = f"{var_theta} = arctan({var_mu})"
    print(solution)

    print("\nIn this equation:")
    print(f"- '{var_theta}' is the angle between the rod and the table surface.")
    print(f"- '{var_mu}' is the coefficient of static friction between the rod and the corner of the table.")

# Execute the function to print the result.
print_sliding_angle_equation()
