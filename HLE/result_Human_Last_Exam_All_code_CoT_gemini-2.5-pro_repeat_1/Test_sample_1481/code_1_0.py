import sympy as sp

def analyze_billiard_generating_function():
    """
    Performs an asymptotic analysis of the billiard generating function H(s, s')
    using symbolic mathematics.
    """
    # Define symbols for the calculation.
    # ds represents the arc-length separation, Δs = s' - s.
    # k represents the curvature at s, κ(s).
    ds, k = sp.symbols('Δs k')

    # Step 1 & 2: Taylor expand the vector q(s') - q(s) in the Frenet frame (T, N) at s.
    # The vector is (q(s') - q(s)) = q'(s)Δs + (1/2)q''(s)Δs² + (1/6)q'''(s)Δs³ + O(Δs⁴).
    # Using Frenet-Serret formulas (q'=T, q''=kN, q'''=-k²T + k'N), we find the
    # components of this vector along the tangent T and normal N.
    # We only need terms that will contribute to H up to O(Δs⁴), so the k' term can be ignored.
    
    # Component along the tangent vector T:
    component_T = ds - k**2 * ds**3 / 6

    # Component along the normal vector N:
    component_N = k * ds**2 / 2

    # Step 3: Compute the squared distance H(s, s')² = |q(s') - q(s)|²
    H_squared = component_T**2 + component_N**2
    
    # Expand the squared expression up to the O(Δs⁵) term to get the correct result for H.
    H_squared_series = sp.series(H_squared, ds, 0, 5).removeO()

    # Step 4: Take the square root and find its series expansion for H(s, s').
    # We assume Δs > 0 for the expansion, but the result holds for |Δs|.
    H = sp.sqrt(H_squared_series)
    H_series = sp.series(H, ds, 0, 5)

    # Format the final expression for clear presentation.
    # The symbol 'ds' is replaced with '|s' - s|' and 'k' with 'κ(s)'.
    final_expr_str = str(H_series.removeO()).replace('Δs', "|s' - s|").replace('k', 'κ(s)')
    
    print("Asymptotic Analysis of Billiard Generating Function H(s, s')")
    print("=" * 60)
    print("The analysis characterizes H(s, s') in the limit |s' - s| -> 0.")
    print(f"\nThe squared distance H(s,s')² expands to: {sp.simplify(H_squared_series)}")
    
    print(f"\nThe resulting asymptotic expansion for H(s,s') is:")
    print(f"H(s, s') = {final_expr_str} + O(|s' - s|⁵)")

    print("\nThis result shows that the leading-order behavior is simply the arc-length separation |s' - s|.")
    print("The next term, -(1/24)κ(s)²|s' - s|³, is the crucial first-order correction, which demonstrates")
    print("the influence of the boundary's local curvature κ(s) on the system's dynamics.")

    # Per the instructions, explicitly output the numbers in the final equation.
    print("\nDecomposition of the final equation's terms:")
    print("H(s, s') = c₁ * |s' - s|¹ + c₂ * κ(s)² * |s' - s|³ + ...")
    
    # First term: 1 * |s'-s|^1
    c1 = H_series.coeff(ds, 1)
    p1 = 1
    print(f"\nTerm 1:")
    print(f"  - Coefficient (c₁): {c1}")
    print(f"  - Power of |s' - s|: {p1}")

    # Second term: -(1/24) * κ(s)^2 * |s'-s|^3
    term2_full_coeff = H_series.coeff(ds, 3) # This is -k**2/24
    c2_numerical = term2_full_coeff.as_coeff_Mul()[0] # This is -1/24
    num, den = c2_numerical.as_numer_denom()
    p2_k = 2
    p3_ds = 3
    print(f"\nTerm 2:")
    print(f"  - Numerical Coefficient (c₂): {num}/{den}")
    print(f"  - Power of curvature κ(s): {p2_k}")
    print(f"  - Power of |s' - s|: {p3_ds}")

analyze_billiard_generating_function()