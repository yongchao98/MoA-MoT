def print_susceptibility_formula():
    """
    This function prints the derived symbolic expression for the magnetic susceptibility (chi)
    for the Ising model on a sparse random graph.
    The expression is built from its symbolic components and then displayed.
    """
    # Define symbolic representations for the physical parameters
    beta = "β"
    c = "c"
    J = "J"

    # Construct the terms of the final equation as strings
    tanh_term = f"tanh({beta}*J)"
    numerator_str = f"{beta} * {c} * {tanh_term}"
    
    denominator_sub_term_str = f"({c} - 1) * {tanh_term}"
    denominator_str = f"1 - {denominator_sub_term_str}"

    # Print the full, final equation for chi
    print("The derived expression for the magnetic susceptibility χ is:")
    print(f"χ = ({numerator_str}) / ({denominator_str})")
    print("\n--- Equation Breakdown ---")

    # "output each number in the final equation" is interpreted as "output each symbolic term"
    print("Numerator terms:")
    print(f"1. Boltzmann factor component: {beta}")
    print(f"2. Graph connectivity: {c}")
    print(f"3. Coupling term: {tanh_term}")
    
    print("\nDenominator terms:")
    print(f"1. Constant: 1")
    print(f"2. Term to be subtracted: {denominator_sub_term_str}")

# Execute the function to display the result
print_susceptibility_formula()