import sympy
from sympy import symbols, Function, series, sqrt, Abs, latex

def analyze_generating_function():
    """
    Performs an asymptotic analysis of the billiard map generating function H(s, s').
    """

    print("Deriving the asymptotic expansion for the billiard generating function H(s, s').\n")
    
    # --- Step 1: Define Symbols ---
    # s is the arc-length parameter
    # ds represents the small separation (s' - s)
    s, ds = symbols('s ds', real=True)
    
    # k(s) is the curvature function kappa(s) at point s
    k = Function('k')(s)
    
    # --- Step 2: Taylor Expansion of the Vector q(s') - q(s) ---
    # We work in the local Frenet-Serret frame (T, N) at point s.
    # In this frame, T = [1, 0] and N = [0, 1].
    #
    # The derivatives of the position vector q are:
    # q'(s)     = T
    # q''(s)    = k(s) * N
    # q'''(s)   = -k(s)**2 * T + k(s).diff(s) * N
    #
    # The vector difference q(s+ds) - q(s) is expanded in a Taylor series:
    # q_diff = q'(s)*ds + (1/2)*q''(s)*ds**2 + (1/6)*q'''(s)*ds**3 + O(ds**4)

    # Component of q_diff along the Tangent vector T
    q_diff_T = ds - (k**2 / 6) * ds**3
    
    # Component of q_diff along the Normal vector N
    q_diff_N = (k / 2) * ds**2

    print("The vector difference q(s') - q(s) in the local (T, N) frame is approximately:")
    print(f"Δq ≈ [ {latex(q_diff_T)} ] T + [ {latex(q_diff_N)} ] N + O(Δs⁴)\n")

    # --- Step 3: Calculate the Squared Norm ||q(s') - q(s)||² ---
    # The squared norm is the sum of the squares of the components.
    norm_sq = q_diff_T**2 + q_diff_N**2
    
    # We expand this expression as a series in ds around ds=0 up to O(ds^5)
    # The expansion needs to go to ds^4 to get the first correction term for the distance.
    norm_sq_series = series(norm_sq, ds, 0, 5)

    print("The squared distance ||q(s') - q(s)||² is expanded as:")
    print(f"||Δq||² ≈ {latex(norm_sq_series)}\n")

    # --- Step 4: Calculate the Norm ||q(s') - q(s)|| ---
    # The distance is the square root of the squared norm.
    # We use Abs(ds) because sqrt(ds**2) = |ds|.
    
    # To correctly expand the square root, factor out ds**2
    norm_sq_poly = norm_sq_series.removeO()
    term_inside_sqrt = (norm_sq_poly / ds**2).expand()
    
    # The distance is |ds| * sqrt(...)
    distance = Abs(ds) * sqrt(term_inside_sqrt)
    
    # Now, find the series expansion for the distance
    distance_series = series(distance, ds, 0, 5)

    print("Taking the square root, the distance ||q(s') - q(s)|| is:")
    print(f"||Δq|| ≈ {latex(distance_series)}\n")
    
    # --- Step 5: Final Asymptotic Form of H(s, s') ---
    # H(s, s') is the negative of the distance.
    H_series = -distance_series
    H_poly = H_series.removeO()

    print("Finally, the generating function H(s, s') = -||Δq|| has the asymptotic form:")
    # Replace ds with |s' - s| for a clearer physical representation
    # Note: ds**2 is Abs(ds)**2, so we can replace ds**2 with |s' - s|**2
    final_expr_str = str(H_poly).replace("ds", "|s' - s|").replace("Abs(|s' - s|)", "|s' - s|")
    
    # To format it in a more standard way for polynomials of |s'-s|
    k_sq = symbols('k(s)**2')
    s_sep = symbols("|s'-s|")
    final_expr = -s_sep + k_sq/24 * s_sep**3
    
    print(f"H(s, s') ≈ {latex(final_expr)}\n")

    # --- Step 6: Output Coefficients as Requested ---
    print("Characterizing the leading-order behavior:")
    print("-" * 40)
    c1_num, c1_den = H_poly.coeff(Abs(ds), 1).as_numer_denom()
    c2_num, c2_den = H_poly.coeff(Abs(ds)**3, 1).as_numer_denom()

    print(f"The equation is of the form: C₁ * |s' - s| + C₂ * k(s)² * |s' - s|³")
    print(f"The coefficient for the leading term, C₁, is: {c1_num}")
    print(f"The coefficient for the first curvature-dependent correction term, C₂, is: {c2_num}/{c2_den}")
    print("-" * 40)

if __name__ == "__main__":
    analyze_generating_function()
