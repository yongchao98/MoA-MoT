import sympy

def analyze_curve_reduction(label, poly_expr):
    """
    Analyzes the reduction property of a curve y^2 = f(x) at a small prime p > 2.

    For a curve defined by the polynomial f(x), this function first finds the primes of
    bad reduction by factoring the discriminant of f(x). It then chooses a small prime 'p'
    of good reduction.

    Finally, it determines if the reduction at 'p' is ordinary or supersingular by
    calculating the Hasse-Witt determinant. For the genus 2 curves in this problem,
    the reduction is ordinary if this determinant is non-zero mod p.
    """
    x = sympy.Symbol('x')
    
    # The equation for the curve is z^2 = f(x)
    print(f"Curve {label}: z^2 = {poly_expr}")

    # Calculate the discriminant to find primes of bad reduction
    try:
        disc = sympy.discriminant(poly_expr, x)
        if disc == 0:
            bad_primes = "N/A (Polynomial has repeated roots)"
        else:
            bad_primes = sorted([p for p in sympy.factorint(disc).keys() if p > 2])
    except Exception as e:
        bad_primes = f"Could not be determined ({e})"

    print(f"  - Primes of bad reduction (> 2): {bad_primes}")

    # Find a small prime p > 2 for testing (must be a prime of good reduction)
    test_p = 3
    if isinstance(bad_primes, list):
        while test_p in bad_primes:
            test_p = sympy.nextprime(test_p)
    
    print(f"  - Testing for reduction type at prime p = {test_p}.")

    # The genus g for polynomials of degree 5 or 6 is 2.
    # The Hasse-Witt test for genus 2 is used.
    # We work in the finite field F_p.
    poly_fp = sympy.Poly(poly_expr, x, modulus=test_p)

    # We need to compute f(x)^((p-1)/2) mod p.
    power_poly = poly_fp ** ((test_p - 1) // 2)

    # The Hasse-Witt determinant for g=2 is: c_{p-1}*c_{2p-2} - c_{p-2}*c_{2p-1}
    # where c_i is the coefficient of x^i in the powered-up polynomial.
    c = lambda i: power_poly.coeff_monomial(x**i)
    det_W = (c(test_p - 1) * c(2*test_p - 2) - c(test_p - 2) * c(2*test_p - 1)) % test_p
    
    # An ordinary reduction corresponds to a non-zero determinant.
    reduction_type = "Ordinary" if det_W != 0 else "Supersingular"
    
    print(f"  - The Hasse-Witt determinant mod {test_p} is {det_W}.")
    print(f"  - Conclusion: The curve has {reduction_type} reduction at p = {test_p}.\n")

if __name__ == '__main__':
    # Define the variable and the polynomials for each curve from the answer choices
    x = sympy.Symbol('x')
    curves = {
        "A": x**5 + 3,
        "B": x**5 - 1,
        "C": x**6 - 1,
        "D": 2*x**5 + 2*x**3 + 1,
        "E": 4*x**5 + 4*x**3 + x**2 + 4*x,
    }

    # Analyze each curve
    for label, poly in curves.items():
        analyze_curve_reduction(label, poly)
