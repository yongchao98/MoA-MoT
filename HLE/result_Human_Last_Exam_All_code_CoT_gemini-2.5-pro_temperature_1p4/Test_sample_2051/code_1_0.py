def calculate_volume_polynomial_degree():
    """
    Calculates the degree of the volume polynomial Z_g,n.

    The degree of the Weil-Petersson volume polynomial Z_g,n as a function
    of boundary lengths is given by the formula: 6g - 6 + 2n.
    """
    # Parameters from the question
    g = 0      # Genus
    n_plus = 3 # Number of positively oriented boundaries
    n_minus = 1 # Number of negatively oriented boundaries

    # Total number of boundaries
    n = n_plus + n_minus

    # Calculate the degree using the formula
    degree = 6 * g - 6 + 2 * n

    # Output the explanation and the final equation
    print(f"For g = {g}, n+ = {n_plus}, and n- = {n_minus}:")
    print(f"The total number of boundaries is n = n+ + n- = {n_plus} + {n_minus} = {n}.")
    print("The degree of the polynomial Z_g,n is given by the formula: 6g - 6 + 2n.")
    print("Substituting the values:")
    print(f"Degree = 6 * {g} - 6 + 2 * {n} = {6 * g} - 6 + {2 * n} = {degree}")
    print("\n---\n")
    print("Final Answer to the question:")
    print("(a) The property of being the volume of a moduli space ensures the function is continuous. So, Yes.")
    print(f"(b) The degree of the polynomial is {degree}.")


calculate_volume_polynomial_degree()