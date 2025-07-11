import sys

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a given Calabi-Yau Link.

    The invariant c is defined as c = d - Σw_i, where d is the degree of the
    quasi-homogeneous polynomial and w_i are the weights of the variables.
    """
    # Weights for the variables (z1, z2, z3, z4, z5)
    weights = [22, 29, 49, 50, 75]

    # The polynomial is:
    # 0 = z1^8*z3 + z1^4*z2^3*z3 + z1*z2^7 + z1*z2*z3*z4*z5 + z2*z3^4 + z4^3*z5 + z5^3
    # For a Calabi-Yau link defined by a quasi-homogeneous polynomial, all
    # monomials must have the same weighted degree.
    # Let's calculate the degree 'd' from a consistent term, e.g., z1^8*z3.
    # Note: The term z1^4*z2^3*z3 has degree 224, which is likely a typo.
    # All other terms have degree 225. We proceed with d=225.
    d = 8 * weights[0] + 1 * weights[2]

    # Calculate the sum of the weights
    sum_of_weights = sum(weights)

    # Calculate the Crawley-Nordström invariant
    c = d - sum_of_weights

    # Print the explanation and the final equation
    print(f"The Crawley-Nordström invariant `c` is calculated using the formula: c = d - Σw_i")
    print(f"Given weights: w = {tuple(weights)}")
    print(f"The quasi-homogeneous degree is d = {d}.")
    print(f"The sum of weights is Σw_i = {' + '.join(map(str, weights))} = {sum_of_weights}.")
    print("\nThe final calculation is:")
    print(f"c = {d} - ({' + '.join(map(str, weights))})")
    print(f"c = {d} - {sum_of_weights}")
    print(f"c = {c}")

    # Output the final answer in the specified format for parsing
    # This is a special marker for the platform, not for the user to copy.
    # The actual numerical answer is printed above for the user.
    sys.stdout.write(f"\n<<<{c}>>>")

if __name__ == '__main__':
    calculate_crawley_nordstrom_invariant()