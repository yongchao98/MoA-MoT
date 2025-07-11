import math

def print_asymptotic_expansion():
    """
    This function presents the asymptotic analysis of the billiard generating function H(s, s').
    It programmatically constructs and prints the formula for the leading-order behavior
    in the limit as the arc-length separation |s' - s| approaches zero.
    """

    # Define the symbolic components of the equation for clarity
    generating_function = "H(s, s')"
    separation = "|s' - s|"
    curvature = "κ(s)"

    # Define the numerical constants from the derivation
    term1_coefficient = 1
    term2_numerator = 1
    term2_denominator = 24
    
    # Define the exponents from the derivation
    sep_exp_term1 = 1
    curv_exp_term2 = 2
    sep_exp_term2 = 3
    big_o_exp = 5

    # --- Explanation Header ---
    print("This script provides the result of the asymptotic analysis for the planar billiard generating function H(s,s').")
    print("-" * 80)
    print("The analysis characterizes the behavior of H(s,s') as the collision points s and s' on the")
    print("billiard boundary become infinitesimally close. H(s,s') represents the Euclidean distance (chord length)")
    print("between these points.")
    print("\nThe derived leading-order expansion, incorporating the boundary's local curvature κ(s), is:\n")

    # --- Print the Final Equation ---
    # The f-string constructs the equation from the predefined components.
    # This fulfills the requirement to output each number in the final equation.
    print(
        f"{generating_function} = "
        f"({term1_coefficient} * {separation}^{sep_exp_term1}) - "
        f"({term2_numerator}/{term2_denominator}) * {curvature}^{curv_exp_term2} * {separation}^{sep_exp_term2} + "
        f"O({separation}^{big_o_exp})"
    )

    # --- Interpretation of the Result ---
    print("\n" + "-" * 80)
    print("Interpretation of the terms:")
    print(f"1. Leading Term ({separation}): This is simply the arc-length distance between the two points.")
    print("   In the absence of curvature (a straight line boundary), this would be the exact distance.")
    print(f"2. Correction Term (-({term2_numerator}/{term2_denominator}) * {curvature}² * {separation}³): This is the first geometric correction.")
    print("   - Its negative sign signifies that for any curved boundary (κ ≠ 0), the chord length H(s,s') is always less than the arc length |s'-s|.")
    print(f"   - The dependence on κ(s)² shows that the magnitude of the correction is the same for convex (κ>0) and concave (κ<0) boundaries.")
    print(f"   - The cubic dependence on the separation ({separation}³) demonstrates that curvature effects are negligible for very small separations but become rapidly important as the separation increases.")


# Execute the function to display the analysis and result
print_asymptotic_expansion()

# The final answer is the formula derived and displayed by the script.
# H(s, s') = |s' - s| - (1/24)κ(s)²|s' - s|³ + O(|s' - s|⁵)
final_formula_string = "H(s, s') = |s' - s| - (1/24)κ(s)²|s' - s|³ + O(|s' - s|⁵)"