def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau link.

    The Crawley-Nordström invariant is defined as d - sum(w_i), where d is the
    degree of the polynomial and w_i are the weights of the ambient space.
    For a Calabi-Yau hypersurface, the degree d must equal the sum of the weights.
    """
    # The weights of the ambient space
    weights = [22, 29, 49, 50, 75]

    # For a Calabi-Yau hypersurface, the degree d equals the sum of the weights.
    sum_of_weights = sum(weights)
    d = sum_of_weights

    # Calculate the Crawley-Nordström invariant
    invariant = d - sum_of_weights

    # Format the weights as a string for printing the equation
    weights_str = " + ".join(map(str, weights))

    # Print the final equation with all the numbers
    print(f"The Crawley-Nordström invariant is calculated as d - (w1 + w2 + w3 + w4 + w5).")
    print(f"Given the Calabi-Yau condition, d = sum of weights.")
    print(f"d = {weights_str} = {d}")
    print(f"Final Equation: {d} - ({weights_str}) = {invariant}")

calculate_crawley_nordstrom_invariant()