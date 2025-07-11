def print_force_equation():
    """
    This function prints the derived formula for the force per unit area on the conductor.
    The formula is constructed from its constituent symbolic parts.
    """
    # Define the symbols used in the equation as string variables for printing
    mu_0 = "μ₀"
    mu = "μ"
    K_0 = "K₀"
    a = "a"
    d = "d"
    y = "y"
    i_x = "î_x"

    # Construct the numerator and denominator of the fraction
    numerator = f"{mu_0}/2 * {K_0}² * sin²({a}{y})"
    denominator_base = f"cosh({a}{d}) + ({mu_0}/{mu}) * sinh({a}{d})"
    denominator = f"({denominator_base})²"

    # Print the final equation for the force per unit area
    print("The force per unit y-z area on the x = d interface is:")
    print(f"f/Area = ( {numerator} ) / ( {denominator} ) * {i_x}")

# Execute the function to display the result
print_force_equation()