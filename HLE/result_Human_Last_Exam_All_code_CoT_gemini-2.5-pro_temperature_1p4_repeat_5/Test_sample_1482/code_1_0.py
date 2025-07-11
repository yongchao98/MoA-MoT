def display_it3_formulation():
    """
    This script prints the mathematical formulation for a vertical cross-section
    of an Interval Type-3 (IT3) Gaussian Upper Membership Function (UMF).
    """

    # Define the components of the equation string for clarity
    base_gaussian = "exp( -0.5 * ( (x' - c) / σ_bar )^2 )"
    type3_exponent = "α_upper^(1 - u)"

    # Construct the full equation string
    full_equation = f"μ_upper(x', u) = [ {base_gaussian} ] ^ ( {type3_exponent} )"

    print("### Mathematical Formulation for a Vertical Cross-Section of an IT3 UMF ###")
    print("\n" + "="*80)
    print("The general mathematical formulation is:")
    print(f"\n    {full_equation}\n")
    print("="*80)

    print("\nWhere the components of the equation are:")
    print("-" * 40)

    # Print each component and its explanation
    print("μ_upper(x', u):")
    print("    The membership value of the vertical cross-section for the Upper")
    print("    Membership Function (UMF).\n")

    print("x':")
    print("    A fixed, constant value of the primary input variable x.\n")

    print("u:")
    print("    The secondary input variable, where u is in the interval [0, 1].")
    print("    It parameterizes the 'blurriness' of the uncertainty bounds.\n")

    print("c:")
    print("    The center (or mean) of the base Gaussian function.\n")

    print("σ_bar (sigma_bar):")
    print("    The standard deviation for the upper bound of the principal Type-2")
    print("    membership function. It defines the widest spread of the Gaussian curve.\n")

    print("α_upper (alpha_upper):")
    print("    A parameter in the interval [0, 1] that controls the degree of blur")
    print("    for the UMF. A value closer to 0 indicates more blurriness.\n")

    print("exp(...):")
    print("    The exponential function (e raised to the power of the argument).\n")

if __name__ == "__main__":
    display_it3_formulation()
