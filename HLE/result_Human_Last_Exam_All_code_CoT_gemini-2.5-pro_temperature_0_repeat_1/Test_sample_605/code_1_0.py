import sys

def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for a given Calabi-Yau link.
    """
    # The weights of the ambient space
    weights = [22, 29, 49, 50, 75]

    # The degree 'd' of the quasi-homogeneous polynomial.
    # For a Calabi-Yau hypersurface, the degree is equal to the sum of the weights.
    # We can also calculate it from the monomials, e.g., for z_5^3:
    # d = 3 * 75 = 225
    # Or for z_1^8*z_3:
    # d = 8 * 22 + 1 * 49 = 176 + 49 = 225
    d = 225

    # Calculate the sum of the weights
    sum_of_weights = sum(weights)

    # The Crawley-Nordström invariant is defined as c = d - sum(weights)
    invariant = d - sum_of_weights

    # Print the components of the calculation
    print(f"The degree of the polynomial is d = {d}.")
    print(f"The sum of the weights is {weights[0]} + {weights[1]} + {weights[2]} + {weights[3]} + {weights[4]} = {sum_of_weights}.")
    
    # Print the final equation and result
    print("\nThe Crawley-Nordström invariant is c = d - \u03A3w\u1D62")
    print(f"c = {d} - {sum_of_weights}")
    print(f"c = {invariant}")
    
    # Output the final answer in the specified format
    sys.stdout.write(f"\n<<<{invariant}>>>\n")

solve_crawley_nordstrom()