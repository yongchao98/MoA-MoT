def solve_moduli_volume_properties():
    """
    Solves the two-part question about the properties of the moduli space volume Z.
    """
    # Part (a): Continuity
    # The function Z_{g, n_+, n_-} represents the volume of a moduli space.
    # While a generic piecewise polynomial function is not guaranteed to be continuous,
    # volumes of well-behaved geometric spaces are continuous functions of their
    # defining parameters. A small perturbation of boundary lengths results in a
    # small perturbation of the volume. Therefore, the function is continuous.
    answer_a = "Yes"
    print(f"(a) {answer_a}")
    print("-" * 20)

    # Part (b): Degree of the polynomial
    # The degree of the polynomial Z_{g, n} is given by the formula: 3g - 3 + n.
    # We are given the case for g = 0, n_+ = 3, and n_- = 1.

    # Define the parameters
    g = 0
    n_plus = 3
    n_minus = 1

    # Calculate the total number of boundaries, n
    n = n_plus + n_minus

    print("(b) To determine the degree of the polynomial Z_{0,3,1}, we use the formula:")
    print("    Degree = 3g - 3 + n")
    print(f"\nGiven values are: g = {g}, n_+ = {n_plus}, n_- = {n_minus}")
    
    # Show calculation for n
    print(f"First, we find the total number of boundaries n:")
    print(f"    n = n_+ + n_- = {n_plus} + {n_minus} = {n}")
    
    # Substitute values into the degree formula and show the calculation
    print("\nNow, we substitute g and n into the degree formula and calculate:")
    print(f"    Degree = 3 * {g} - 3 + {n}")
    
    term1 = 3 * g
    print(f"    Degree = {term1} - 3 + {n}")
    
    term2 = term1 - 3
    print(f"    Degree = {term2} + {n}")

    degree = term2 + n
    print(f"    Degree = {degree}")
    
    answer_b = degree
    print(f"\nThe degree of the polynomial is {answer_b}.")

# Execute the function to print the solution
solve_moduli_volume_properties()