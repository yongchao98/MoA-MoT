# This script requires a SageMath environment to execute.
from sage.all import EllipticCurve, DirichletGroup, QQbar, pari

def solve_elliptic_curve_problem():
    """
    Solves the problem by calculating the rank of an elliptic curve and
    the leading coefficients of its twisted L-functions.
    """
    # Step 1: Define the elliptic curve E and find its rank r.
    # The curve y^2 + y = x^3 - x^2 - 10x - 20 has LMFDB label '49.a1'.
    try:
        E = EllipticCurve('49.a1')
    except Exception:
        # Fallback to defining from coefficients if LMFDB is not available.
        E = EllipticCurve([0, -1, 1, -10, -20])

    # The rank r is the Z-rank of the Mordell-Weil group E(Q).
    r = E.rank()

    # Step 2: Define the primitive cubic Dirichlet characters of conductor 7.
    # We use QQbar as the codomain to handle complex numbers (roots of unity).
    G = DirichletGroup(7, QQbar)
    chars = [chi for chi in G if chi.order() == 3 and chi.is_primitive()]

    # To ensure a deterministic order for chi_1 and chi_2, we sort them
    # based on the imaginary part of their value on a generator of (Z/7Z)*, which is 3.
    chi1, chi2 = chars[0], chars[1]
    if chi1(3).imag() < 0:
        chi1, chi2 = chi2, chi1

    # Step 3: Compute the leading coefficients a and b.
    # The twisted L-functions have analytic rank 1, so the leading coefficient
    # is the first derivative at s=1. We use the lfun function from PARI/GP.
    # pari(E).lfun(chi, s, d) computes the d-th derivative of L(E, s, chi).
    a = pari(E).lfun(chi1, 1, 1)
    b = pari(E).lfun(chi2, 1, 1)

    # Step 4: Calculate the sum r + a + b and print the result.
    # As b is the complex conjugate of a, their sum is real.
    total_sum = r + a + b

    # Format the numbers for the final equation output.
    # We display the real and imaginary parts of a and b.
    a_str = f"({a.real():.4f} + {a.imag():.4f}i)"
    b_str = f"({b.real():.4f} - {abs(b.imag()):.4f}i)"
    
    # The final result is the real part of the sum, rounded.
    final_result_rounded = round(total_sum.real(), 4)

    print("The rank of the Mordell-Weil group E(Q) is r.")
    print("The leading coefficients of the twisted L-functions at s=1 are a and b.")
    print("We need to compute r + a + b.")
    print("\nThe final equation with the computed values is:")
    print(f"{r} + {a_str} + {b_str} = {final_result_rounded}")

# Execute the main function.
solve_elliptic_curve_problem()