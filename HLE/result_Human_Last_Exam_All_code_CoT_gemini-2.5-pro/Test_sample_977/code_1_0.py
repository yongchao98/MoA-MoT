def print_potential_formula():
    """
    This script prints the derived mathematical formula for the electric potential Phi(x, y).
    The derivation is based on solving Laplace's equation with the given boundary conditions.
    The final expression is built from its constituent parts for clarity.
    """

    # The denominator of the potential expression is common for both regions.
    # It arises from solving for the unknown coefficients using the interface boundary conditions.
    denominator = "k * (epsilon_2*cosh(k*a)*sinh(k*b) + epsilon_1*sinh(k*a)*cosh(k*b))"

    # The numerator for the potential in the region 0 < y < a. This is the region
    # specifically asked about in the problem.
    numerator_region_2 = "-sigma_0*sinh(k*b)*sinh(k*(y - a))*sin(k*x)"

    # Print the potential for the requested region: 0 <= y <= a
    print("The electric potential Phi(x, y) in the region 0 <= y <= a is:")
    print(f"Phi(x, y) = ({numerator_region_2}) / ({denominator})")

    # For verification against the multiple-choice options, which provide the complete
    # piecewise function, we also state the potential in the other region.
    numerator_region_1 = "sigma_0*sinh(k*a)*sinh(k*(y + b))*sin(k*x)"

    print("\nThe full piecewise potential function provided in the correct answer choice is:")
    print("Phi(x, y) = ")
    print(f"  For 0 < y < a:    ({numerator_region_2}) / ({denominator})")
    print(f"  For -b < y < 0:   ({numerator_region_1}) / ({denominator})")

    print("\nThis derived solution matches answer choice A.")

if __name__ == '__main__':
    print_potential_formula()