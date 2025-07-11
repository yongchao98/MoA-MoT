import sys

def display_generating_function_asymptotics():
    """
    This function provides a rigorous characterization of the asymptotic behavior
    of the generating function H(s, s') in planar Birkhoff billiard dynamics,
    highlighting the role of local boundary curvature κ(s).
    """

    # Using UTF-8 encoding for mathematical symbols
    if sys.stdout.encoding != 'UTF-8':
        sys.stdout.reconfigure(encoding='UTF-8')

    print("--- Asymptotic Analysis of the Billiard Generating Function H(s, s') ---")
    print("\nIn the theory of planar billiards, the generating function H(s, s') connects")
    print("consecutive points of impact on the boundary, s and s'. It is defined as the")
    print("Euclidean distance between these two points.")
    print("\nIn the limit where the separation |s' - s| is very small, we can expand H(s, s')")
    print("to understand how the local geometry of the boundary affects the dynamics.")
    print("\nThe resulting asymptotic expansion is presented below.")

    # Define the numerical coefficients from the derivation
    leading_term_coeff = 1
    curvature_term_coeff_numerator = -1
    curvature_term_coeff_denominator = 24
    leading_term_power = 1
    curvature_term_power = 3
    error_term_power = 4 # O(|s'-s|⁴) from derivation, but O(|s'-s|⁵) is more precise if we keep next order term in the root. 
    # For H, the next term is O((Δs)^5). But in terms of the structure of the prompt asking for leading order it's more clean to say O((Δs)^4) for now.
    # The term we are showing is Δs^3, so the next order is Δs^4. No, my derivation was H ~ |ds| - k^2/24 |ds|^3 + O(|ds|^5) . Let me check the derivation again.
    # H = |ds| (1 - k^2/12 ds^2 - k k'/12 ds^3)^(1/2)
    # H = |ds| (1 - 1/2(k^2/12 ds^2 + k k'/12 ds^3) + O(ds^4))
    # H = |ds| (1 - k^2/24 ds^2 - k k'/24 ds^3 + O(ds^4))
    # H = |ds| - k^2/24 |ds|^3 - k k'/24 ds^4 + O(|ds|^5)
    # So the error term is indeed O(|s'-s|^4) if we truncate here. Okay, O(|s'-s|^4) is correct.
    
    # Mathematical symbols using unicode
    kappa = "\u03BA"
    squared = "\u00B2"
    cubed = "\u00B3"
    fourth_power = "\u2074"
    approx = "\u2248" # "approximately equal to" symbol

    # Construct the final equation string
    equation = (
        f"H(s, s') {approx} |s' - s| - "
        f"({abs(curvature_term_coeff_numerator)}/{curvature_term_coeff_denominator}) "
        f"{kappa}(s){squared} |s' - s|{cubed} + O(|s' - s|{fourth_power})"
    )

    print("\n" + "="*70)
    print("Final Asymptotic Formula:")
    print(f"\n  {equation}\n")
    print("="*70)


    print("\nWhere:")
    print(f"  H(s, s'): The generating function, i.e., distance between q(s) and q(s').")
    print(f"  s, s':    Arc-length parameters on the billiard boundary.")
    print(f"  {kappa}(s):      The local curvature of the boundary at point s.")
    print(f"  |s' - s|: The small arc-length separation between the points.")

    print("\nInterpretation:")
    print("  - The first term, |s' - s|, is the arc length, representing the trivial case of a flat boundary (κ=0).")
    print("  - The second term shows that for a curved boundary (κ ≠ 0), the chord length H is always")
    print("    *less* than the arc length. This deviation is cubic in the separation |s' - s| and quadratic")
    print(f"    in the curvature {kappa}(s). This is the leading-order geometric correction.")

    print("\nExplicit numerical coefficients in the expansion:")
    print(f"  - Coefficient of the leading term |s' - s|: {leading_term_coeff}")
    print(f"  - Coefficient of the curvature term {kappa}(s){squared}|s' - s|{cubed}: {curvature_term_coeff_numerator}/{curvature_term_coeff_denominator}")

if __name__ == '__main__':
    display_generating_function_asymptotics()
