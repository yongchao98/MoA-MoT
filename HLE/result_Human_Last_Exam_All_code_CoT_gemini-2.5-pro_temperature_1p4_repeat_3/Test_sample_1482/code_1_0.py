import math

def display_it3_formulation():
    """
    This function prints the mathematical formulation for the vertical cross-sections
    of a Gaussian-based Interval Type-3 Membership Function (IT3 MF).
    """

    print("In an Interval Type-3 Fuzzy Logic System, a vertical cross-section at a fixed input 'x' is an Interval Type-2 Membership Function defined over a secondary variable 'u'.")
    print("This cross-section is characterized by its upper and lower membership bounds, derived from the overall IT3 MF bounds.")
    print("-" * 80)
    print("The formulation uses a Gaussian paradigm. The parameters are:")
    print("  'x': The primary input variable.")
    print("  'u': The secondary input variable, representing the membership degree in the third dimension.")
    print("  'c_x': The center (mean) of the Gaussian function in the primary dimension.")
    print("  'σ_bar_x' and 'σ_lower_x': The upper and lower bounds for the standard deviation in the 'x' dimension (σ_bar_x > σ_lower_x).")
    print("  'σ_bar_u' and 'σ_lower_u': The upper and lower bounds for the standard deviation in the 'u' dimension (σ_bar_u > σ_lower_u).")
    print("-" * 80)
    print("Mathematical Formulation for the Bounds:")
    print("\n1. Upper Membership Bound (μ_bar):")
    print("   This function defines the upper surface of the IT3 MF and the Upper Membership Function of any vertical cross-section.")
    print("\n   μ_bar(x, u) = exp(-0.5 * ((x - c_x) / σ_bar_x)^2) * exp(-0.5 * ((u - 1) / σ_bar_u)^2)")
    print("\n" + "=" * 80)
    print("\n2. Lower Membership Bound (μ_lower):")
    print("   This function defines the lower surface of the IT3 MF and the Lower Membership Function of any vertical cross-section.")
    print("\n   μ_lower(x, u) = exp(-0.5 * ((x - c_x) / σ_lower_x)^2) * exp(-0.5 * ((u - 1) / σ_lower_u)^2)")
    print("\n" + "-" * 80)

# Execute the function to display the formulation.
display_it3_formulation()
