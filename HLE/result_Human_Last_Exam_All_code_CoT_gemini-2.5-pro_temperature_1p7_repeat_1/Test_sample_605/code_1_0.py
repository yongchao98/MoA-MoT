import math

def calculate_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for the given Calabi-Yau link.
    """
    # The weights of the ambient space P(w1, w2, w3, w4, w5)
    weights = [22, 29, 49, 50, 75]

    # For a Calabi-Yau hypersurface, the degree 'd' of the defining polynomial
    # is equal to the sum of the weights.
    d = sum(weights)

    # The Crawley-Nordström invariant for a Calabi-Yau threefold simplifies to:
    # CN(X) = (1/12) * Σ((w_i^2 + 1) / w_i) for i=1 to 5.
    # The first term of the general formula becomes zero because (Σw_i - d) = 0.

    # We calculate the sum part of the formula.
    sum_of_terms = 0
    for w in weights:
        term = (w**2 + 1) / w
        sum_of_terms += term
        
    # The invariant is the sum divided by 12.
    cn_invariant = sum_of_terms / 12

    # Print the explanation and the final equation with all numbers.
    print("The Crawley-Nordström invariant for a Calabi-Yau hypersurface is calculated using a simplified formula:")
    print("CN(X) = (1/12) * Σ((w_i^2 + 1) / w_i)")
    print("\nFor the weights W = (22, 29, 49, 50, 75), the calculation is as follows:")

    # Construct and print the full equation with numerical values
    equation_str = "CN(X) = 1/12 * ( "
    terms_list = []
    for w in weights:
        terms_list.append(f"({w}^2 + 1)/{w}")
    equation_str += " + ".join(terms_list)
    equation_str += f" ) = {cn_invariant}"
    
    print("\nFinal Equation:")
    print(equation_str)

# Run the calculation
calculate_crawley_nordstrom()