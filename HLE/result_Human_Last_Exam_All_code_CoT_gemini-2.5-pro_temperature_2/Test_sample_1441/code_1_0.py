def calculate_alpha_for_dimensionality(d):
    """
    Calculates the specific heat critical exponent alpha (α)
    using the first-order epsilon expansion for a scalar field theory (N=1).

    The epsilon expansion provides an approximation valid for dimensionalities d
    that are close to the upper critical dimension, d_u = 4.

    Args:
        d (float or int): The spatial dimensionality, must be less than 4.
    """
    # For a scalar field theory, the number of components of the order parameter is N=1.
    N = 1
    # The upper critical dimension for this model is 4.
    upper_critical_dimension = 4

    if d >= upper_critical_dimension:
        print(f"Error: The epsilon expansion is valid for d < {upper_critical_dimension}.")
        print("For d >= 4, mean-field theory applies and α = 0.")
        return

    # Epsilon (ε) is the small parameter for the expansion.
    epsilon = upper_critical_dimension - d

    # The formula for α to first order in ε is: α = (4 - N) / (N + 8) * ε
    numerator = 4 - N
    denominator = N + 8

    # Calculate the final value for alpha
    alpha = (numerator / denominator) * epsilon

    print(f"For spatial dimensionality d = {d}:")
    print(f"The expansion parameter ε = {upper_critical_dimension} - {d} = {epsilon:.4f}")
    print(f"Using the first-order formula for α = [(4 - N) / (N + 8)] * ε:")
    # Output the final equation with all the numbers, as requested.
    print(f"α ≈ [({4} - {N}) / ({N} + {8})] * {epsilon:.4f}")
    print(f"α ≈ [{numerator} / {denominator}] * {epsilon:.4f}")
    print(f"α ≈ {numerator / denominator:.4f} * {epsilon:.4f}")
    print(f"α ≈ {alpha:.4f}")
    print("-" * 30)

def main():
    """
    Main function to demonstrate the calculation for key dimensionalities.
    """
    print("This script calculates the specific heat critical exponent (α) as a quantitative")
    print("function of spatial dimensionality (d) below the upper critical dimension of 4.\n")

    # Demonstrate for d=3 (a common physical dimension)
    calculate_alpha_for_dimensionality(3)

    # Demonstrate for d=2
    calculate_alpha_for_dimensionality(2)

    # Demonstrate for a value very close to the upper critical dimension
    calculate_alpha_for_dimensionality(3.9)


if __name__ == "__main__":
    main()