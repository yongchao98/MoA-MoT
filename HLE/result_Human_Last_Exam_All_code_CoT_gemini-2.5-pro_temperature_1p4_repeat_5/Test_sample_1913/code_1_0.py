# This script should be run in a SageMath environment.
from sage.all import EllipticCurve, DirichletGroup, QQbar, n, real_part

def solve_elliptic_curve_problem():
    """
    Solves the problem by calculating the rank of an elliptic curve and the
    leading coefficients of its twisted L-functions.
    """
    # 1. Define the Elliptic Curve E from its Weierstrass coefficients
    # E: y^2 + y = x^3 - x^2 - 10x - 20
    E = EllipticCurve([0, -1, 1, -10, -20])

    # 2. Compute the rank r of E(Q)
    # For this curve (LMFDB label '49.a3'), the rank is 0.
    r = E.rank()

    # 3. Define the primitive cubic Dirichlet characters of conductor 7
    G = DirichletGroup(7, QQbar)
    # For prime conductor 7, all non-trivial characters are primitive.
    # We select the characters of order 3.
    cubic_chars = [chi for chi in G if chi.order() == 3]
    chi1 = cubic_chars[0]
    chi2 = cubic_chars[1]

    # 4. Compute the leading coefficients a and b of the twisted L-series at s=1.
    # The method .value(1) computes the leading coefficient of the Taylor series at s=1.
    # For these twists, the analytic rank is 1, so this is L'(1).
    a = E.lseries(chi1).value(1)
    b = E.lseries(chi2).value(1)

    # 5. Calculate the final sum r + a + b
    total_sum = r + a + b
    
    # The result should be real since b is the complex conjugate of a.
    # We take the real part and round to four decimal places.
    final_result_rounded = n(real_part(total_sum)).round(4)
    
    # Output the details of the calculation
    print(f"The Z-rank of E(Q) is r = {r}")
    # Using n() to get a numerical approximation for printing
    a_numeric = n(a)
    b_numeric = n(b)
    print(f"The leading coefficient of L(E, s, chi_1) at s=1 is a = {a_numeric}")
    print(f"The leading coefficient of L(E, s, chi_2) at s=1 is b = {b_numeric}")

    # Output the final equation with the computed values
    print("\nThe final equation is:")
    print(f"{r} + ({a_numeric}) + ({b_numeric}) = {final_result_rounded}")

solve_elliptic_curve_problem()