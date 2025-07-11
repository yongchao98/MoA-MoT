def generate_it3_mf_formulation():
    """
    Generates and prints the mathematical formulation for the upper bound
    of a vertical cross-section of a Gaussian Interval Type-3 Membership Function (IT3 MF).
    """

    # The formulation for the upper bound (mu_bar_A) is based on a Gaussian function.
    # In an IT3 MF, a parameter like the standard deviation (sigma) is uncertain and
    # represented by an interval [sigma_lower(u), sigma_bar(u)] that depends on
    # the secondary variable u.

    # The Gaussian function value increases as its standard deviation increases
    # (for any x not equal to the center c).
    # Thus, the upper membership bound (mu_bar) corresponds to the upper
    # standard deviation bound (sigma_bar).
    upper_bound_formula = "mu_bar_A(x, u) = exp(-0.5 * ((x - c) / sigma_bar(u))^2)"

    print("The mathematical formulation that characterizes the upper bound of the vertical cross-section of a Gaussian IT3 MF is:")
    print(upper_bound_formula)
    print("\nWhere the terms are defined as:")
    print(" - mu_bar_A(x, u): The upper membership bound for a primary input 'x' and secondary input 'u'.")
    print(" - exp(...): The exponential function.")
    print(" - x: The primary input variable, which is a fixed value for any given vertical cross-section.")
    print(" - c: The center (or mean) of the Gaussian function.")
    print(" - u: The secondary variable, which typically represents a level of primary membership.")
    print(" - sigma_bar(u): The function defining the upper bound of the standard deviation for a given 'u'.")
    print("\nThe numbers explicitly defined in the final equation are -0.5 and 2 (from the square).")


if __name__ == "__main__":
    generate_it3_mf_formulation()