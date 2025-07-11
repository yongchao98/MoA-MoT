import math

def display_it3_vertical_slice_formulation():
    """
    This function formulates and prints the mathematical equation for a vertical
    cross-section of an Interval Type-3 Membership Function (IT3 MF)
    using a Gaussian-based paradigm.
    """

    # --- Introduction and Definition of Terms ---
    print("In an Interval Type-3 Fuzzy Logic System, the vertical cross-section at a fixed primary input 'x' is an Interval Type-2 Membership Function.")
    print("This function is defined over a secondary domain 'u'.")
    print("The uncertainty in this cross-section is captured by its upper and lower membership functions.")
    print("\nThe mathematical formulation for these bounds is often based on a Gaussian function, where the parameters of the Gaussian depend on 'x'.\n")

    print("--- Symbolic Definitions ---")
    print("Let the formulation for a bound (either upper or lower) of the vertical cross-section be denoted as 'μ_bound(x, u)'.")
    print("  x: The primary input variable.")
    print("  u: The secondary input variable, which lies in the primary membership domain J_x.")
    print("  C(x): The center (mean) of the Gaussian function, which varies as a function of 'x'.")
    print("  Σ(x): The standard deviation (sigma) of the Gaussian function, which also varies as a function of 'x', defining the spread of uncertainty.\n")

    # --- Mathematical Formulation ---
    print("--- Formulation of the Gaussian-based Vertical Cross-Section Bound ---")

    # The components of the equation
    term_1 = "μ_bound(x, u)"
    term_2 = "exp("
    number_1 = -0.5
    term_3 = "* ( (u - C(x)) / Σ(x) )"
    number_2 = 2
    term_4 = ")"

    # Print the final equation piece by piece to highlight each number
    print("The equation is constructed as follows:")
    print(f"  {term_1} = {term_2}{number_1} {term_3}^{number_2} {term_4}")
    print("\nWhere:")
    print(f"  - The exponential function 'exp' defines the Gaussian shape.")
    print(f"  - The constant '{number_1}' and the exponent '{number_2}' are standard components of a Gaussian probability density function.")
    print(f"  - 'u - C(x)' calculates the distance of the secondary variable from the center of the distribution for a given 'x'.")
    print(f"  - 'Σ(x)' normalizes this distance by the standard deviation for that given 'x'.")


if __name__ == '__main__':
    display_it3_vertical_slice_formulation()

<<<μ_bound(x, u) = exp(-0.5 * ( (u - C(x)) / Σ(x) )^2)>>>