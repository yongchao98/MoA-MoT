import sympy as sp

def billiard_generating_function_asymptotics():
    """
    Performs and explains the asymptotic analysis of the billiard generating function H(s, s').

    This function uses the sympy library to derive the leading-order behavior of H(s,s')
    in the limit as the arc-length separation |s'−s| approaches zero, including the
    influence of the boundary's local curvature κ(s).
    """

    # --- Introduction and Setup ---
    print("Deriving the asymptotic expansion for the billiard generating function H(s, s')")
    print("=" * 70)
    print("The generating function is the Euclidean distance between boundary points s and s'.")
    print("H(s, s') = |r(s') - r(s)|\n")

    # Define symbolic variables
    # ds represents the small separation, Δs = s' - s
    # kappa_s represents the local curvature κ(s) at point s
    ds = sp.Symbol("Δs")
    kappa_s = sp.Symbol("κ(s)")

    # --- Step 1: Taylor Expansion of the Displacement Vector ---
    print("Step 1: Define the displacement vector r(s') - r(s) using a Taylor series in a local Frenet-Serret frame (T, N).")
    print("In this frame, r(s) is the origin, T=<1,0>, and N=<0,1>.")
    print("r(s') - r(s) ≈ r'(s)Δs + (1/2!)r''(s)Δs² + (1/6)r'''(s)Δs³")
    print(f"where r'(s)=T, r''(s)=κ(s)N, and r'''(s)=-κ(s)²T + κ'(s)N.\n")

    # The κ'(s) term's contribution to H is of order O(Δs⁵), so it can be neglected for this analysis.
    # Displacement vector dr = r(s') - r(s) components
    dr_x = 1 * ds - sp.Rational(1, 6) * kappa_s**2 * ds**3
    dr_y = sp.Rational(1, 2) * kappa_s * ds**2

    print(f"Displacement x-component (along T): {dr_x}")
    print(f"Displacement y-component (along N): {dr_y}\n")

    # --- Step 2: Calculate the Squared Distance H² ---
    print("Step 2: Calculate the squared distance H² = |r(s') - r(s)|².")
    H_squared = sp.expand(dr_x**2 + dr_y**2)

    # Use series expansion to keep terms up to a certain order, simplifying the view
    H_squared_series = H_squared.series(ds, 0, 5)

    print(f"H² = (x-component)² + (y-component)²")
    print(f"   ≈ {H_squared_series.removeO()}\n")

    # --- Step 3: Find H by taking the Square Root ---
    print("Step 3: Calculate H = sqrt(H²) and find its series expansion for small Δs.")
    H = sp.sqrt(H_squared)
    H_series = H.series(ds, 0, 5)

    print(f"H ≈ {H_series}\n")

    # --- Final Result ---
    print("Step 4: Final Asymptotic Expression.")
    print("=" * 70)
    print("The asymptotic expansion of the generating function H(s, s'), replacing Δs with |s'-s|, is:")

    final_expression = H_series.removeO()
    c1 = final_expression.coeff(ds, 1)
    c3_raw = final_expression.coeff(ds, 3)

    # Use sympy.pretty to format the output nicely
    eq_str = sp.pretty(c1 * sp.Symbol("|s'-s|") + c3_raw * sp.Symbol("|s'-s|")**3, use_unicode=False)
    final_eq = f"H(s, s') = {eq_str} + O(|s'-s|⁵)"
    print(final_eq)

    print("\nThis shows the leading behavior is the arc length, with a third-order correction from curvature.")

    print("\nThe coefficients for each term in the equation are:")
    c3_coeff = c3_raw / kappa_s**2
    print(f"Term 1 (|s'-s|): Coefficient = {c1}")
    print(f"Term 2 (κ(s)²|s'-s|³): Coefficient = {c3_coeff}")
    print("=" * 70)


if __name__ == "__main__":
    billiard_generating_function_asymptotics()
<<<H(s, s') = |s' - s| - κ(s)²⋅|s' - s|³/24 + O(|s' - s|⁵)>>>