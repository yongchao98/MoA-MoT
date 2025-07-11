def print_chromatic_polynomial():
    """
    Prints the chromatic polynomial for the given graph.
    The polynomial P(k) = k^5 - 7k^4 + 19k^3 - 23k^2 + 10k was derived
    using the deletion-contraction method or direct counting.
    This function formats and prints the final equation, including all numbers.
    """
    
    # The chromatic polynomial equation as a string.
    # Each number (coefficients and exponents) is explicitly part of the string.
    poly_string = "P(k) = k^5 - 7*k^4 + 19*k^3 - 23*k^2 + 10*k"
    
    print("The chromatic polynomial for the graph is:")
    print(poly_string)

print_chromatic_polynomial()