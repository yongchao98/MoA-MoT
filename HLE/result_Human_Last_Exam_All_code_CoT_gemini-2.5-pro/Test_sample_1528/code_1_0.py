# This script is intended to be run in a SageMath environment.
# You can use the free online SageMathCell: https://sagecell.sagemath.org/

try:
    # Import necessary functions from the SageMath library
    from sage.all import PolynomialRing, QQ, HyperellipticCurve, factor
except ImportError:
    print("This script must be run in a SageMath environment.")
else:
    # The curve is defined by the equation y^2 = f(x).
    # We first define the polynomial f(x) over the field of rational numbers (QQ).
    R, x = PolynomialRing(QQ, 'x').objgen()
    # The polynomial is y^2 = 1*x^6 + 4*x^5 + 6*x^4 + 2*x^3 + 1*x^2 + 2*x + 1
    f = x**6 + 4*x**5 + 6*x**4 + 2*x**3 + x**2 + 2*x + 1

    # For clarity, let's print the equation's coefficients.
    print("Curve equation: y^2 = (1)*x^6 + (4)*x^5 + (6)*x^4 + (2)*x^3 + (1)*x^2 + (2)*x + (1)")

    # Create the Hyperelliptic Curve object from the polynomial.
    C = HyperellipticCurve(f)

    # The conductor of a hyperelliptic curve is the conductor of its Jacobian variety.
    J = C.jacobian()

    # Compute the conductor.
    N = J.conductor()

    # Print the final result.
    print(f"\nThe conductor of the curve is: {N}")

    # To satisfy the requirement of outputting numbers in the final equation,
    # we show the prime factorization of the conductor.
    N_factors = factor(N)

    # Format the factorization into a human-readable string: N = p1 * p2 * ...
    factor_list = [str(p) for p, e in N_factors for _ in range(e)]
    factor_equation = " * ".join(factor_list)

    print(f"The prime factorization of the conductor is {N} = {factor_equation}")
