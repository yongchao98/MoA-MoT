def calculate_polynomial_degree():
    """
    Calculates the degree of the polynomial Z_{g, n_+, n_-} for the specified case.
    """
    # Parameters for the specific case
    g = 0
    n_plus = 3
    n_minus = 1

    # The degree of the polynomial Z_{g, n_+, n_-} is given by the formula:
    # degree = 3g - 3 + n_+ + n_-
    
    # Perform the calculation
    degree = 3 * g - 3 + n_plus + n_minus

    print("Calculation for the degree of the polynomial Z_{g, n_+, n_-}:")
    print(f"Given parameters: g = {g}, n_+ = {n_plus}, n_- = {n_minus}")
    print("The degree is calculated using the formula: 3*g - 3 + n_+ + n_-")
    print("Substituting the values into the formula:")
    print(f"Degree = 3*({g}) - 3 + {n_plus} + {n_minus} = {degree}")

calculate_polynomial_degree()