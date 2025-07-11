import sys
from collections import Counter

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordstrom invariant for a given Calabi-Yau link.
    """
    # Weights for the variables (z1, z2, z3, z4, z5)
    weights = [22, 29, 49, 50, 75]

    # The polynomial is f = z1^8*z3 + z1^4*z2^3*z3 + z1*z2^7 + z1*z2*z3*z4*z5 + z2*z3^4 + z4^3*z5 + z5^3
    # We represent each term by the tuple of its exponents.
    exponents = [
        (8, 0, 1, 0, 0),  # z1^8*z3
        (4, 3, 1, 0, 0),  # z1^4*z2^3*z3
        (1, 7, 0, 0, 0),  # z1*z2^7
        (1, 1, 1, 1, 1),  # z1*z2*z3*z4*z5
        (0, 1, 4, 0, 0),  # z2*z3^4
        (0, 0, 0, 3, 1),  # z4^3*z5
        (0, 0, 0, 0, 3)   # z5^3
    ]

    # Calculate the weighted degree for each term
    term_degrees = []
    for exp in exponents:
        degree = sum(e * w for e, w in zip(exp, weights))
        term_degrees.append(degree)

    # For a weighted homogeneous polynomial, all term degrees must be equal.
    # We check this and find the polynomial's degree 'd'.
    degree_counts = Counter(term_degrees)
    if len(degree_counts) > 1:
        # This indicates a potential typo in the polynomial, as most terms
        # for a Calabi-Yau hypersurface have the same weighted degree.
        # We proceed with the most frequent degree.
        print(f"Warning: The polynomial is not perfectly weighted-homogeneous. The degrees of its terms are {term_degrees}.", file=sys.stderr)
        d = degree_counts.most_common(1)[0][0]
        print(f"Using the most common degree, d = {d}\n")
    else:
        d = term_degrees[0]
        print(f"The polynomial is weighted-homogeneous with degree d = {d}\n")

    # Calculate the sum of the weights
    sum_of_weights = sum(weights)

    # Calculate the Crawley-Nordstrom invariant
    # Formula: c(f) = (d - sum_of_weights) / d
    numerator = d - sum_of_weights
    denominator = d

    if denominator == 0:
      print("Error: The degree 'd' is zero, so the invariant cannot be calculated.", file=sys.stderr)
      return

    invariant = numerator / denominator
    # If the result is a whole number, display it as an integer.
    if invariant == int(invariant):
        invariant = int(invariant)
    
    # Print the final calculation as requested
    print("The Crawley-Nordstrom invariant is calculated using the formula:")
    print("c(f) = (d - Σwᵢ) / d")
    print("\nStep 1: The degree 'd' is the weighted degree of the polynomial.")
    print(f"d = {d}")
    print("\nStep 2: The sum of weights Σwᵢ is:")
    print(f"Σwᵢ = {' + '.join(map(str, weights))} = {sum_of_weights}")
    print("\nStep 3: Plugging the numbers into the formula:")
    print(f"c(f) = ({d} - {sum_of_weights}) / {d}")
    print(f"c(f) = {numerator} / {denominator}")
    print("\nFinal Result:")
    print(f"The Crawley-Nordstrom invariant is: {invariant}")

if __name__ == '__main__':
    calculate_crawley_nordstrom_invariant()