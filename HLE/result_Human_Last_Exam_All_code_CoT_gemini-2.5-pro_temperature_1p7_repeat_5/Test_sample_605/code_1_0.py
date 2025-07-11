import re

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a given Calabi-Yau link.

    The invariant is the sum of the weights of variables that appear as
    pure power monomials (e.g., z_i^k) in the defining polynomial.
    """
    # The polynomial defining the Calabi-Yau Link
    polynomial = "z_1^8z_3+z_1^4z_2^3z_3+z_1z_2^7+z_1z_2z_3z_4z_5+z_2z_3^4+z_4^3z_5+z_5^3"

    # The weights of the variables (z1, z2, z3, z4, z5)
    weights = {
        '1': 22,
        '2': 29,
        '3': 49,
        '4': 50,
        '5': 75
    }

    # Split the polynomial into its monomial terms
    terms = polynomial.split('+')

    # A regular expression to find all variable indices (the 'i' in z_i) in a term
    variable_pattern = re.compile(r'z_(\d+)')

    contributing_weights = []
    pure_power_variables = []

    # Analyze each term to see if it's a pure power of a single variable
    for term in terms:
        # Find all variables present in the current term
        variables_in_term = set(variable_pattern.findall(term))

        # A pure power monomial will only contain one unique variable
        if len(variables_in_term) == 1:
            var_index = variables_in_term.pop()
            contributing_weights.append(weights[var_index])
            pure_power_variables.append(f"z_{var_index}")
    
    # Calculate the invariant by summing the collected weights
    invariant_value = sum(contributing_weights)

    # Print the explanation and the final result in an equation format
    print("The Crawley-Nordström invariant is the sum of the weights of variables")
    print("that appear as pure power monomials in the polynomial.")
    print("-" * 50)
    if not contributing_weights:
        print("No pure power monomials found.")
        print("The invariant is 0.")
        print("0 = 0")
    else:
        weights_str = " + ".join(map(str, contributing_weights))
        variables_str = ", ".join(pure_power_variables)
        print(f"The pure power monomial(s) correspond to the variable(s): {variables_str}.")
        print(f"The respective weight(s) are: {', '.join(map(str, contributing_weights))}.")
        print("\nThe resulting sum is:")
        print(f"{weights_str} = {invariant_value}")

# Run the calculation
calculate_crawley_nordstrom_invariant()