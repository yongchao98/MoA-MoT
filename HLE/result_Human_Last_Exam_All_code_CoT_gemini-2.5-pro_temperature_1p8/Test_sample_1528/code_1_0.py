# This is a SageMath script. It must be run in a SageMath environment.

def find_conductor():
    """
    This function computes the conductor of a given hyperelliptic curve.
    """
    # Define the polynomial ring over the rational numbers QQ
    R_x = PolynomialRing(QQ, 'x')
    x = R_x.gen()

    # Define the polynomial f(x) from the equation y^2 = f(x)
    # The coefficients are 1, 4, 6, 2, 1, 2, 1
    f = x**6 + 4*x**5 + 6*x**4 + 2*x**3 + x**2 + 2*x + 1

    # Create the hyperelliptic curve object
    C = HyperellipticCurve(f)

    # The conductor is an attribute of the Jacobian of the curve.
    # We get the Jacobian of C.
    J = C.jacobian()

    # Compute the conductor. SageMath has a built-in function for this.
    N = J.conductor()

    # The final instruction asks to output each number in the final equation.
    # We will print a sentence containing the original equation and the computed conductor.
    print(f"The conductor of the curve defined by y^2 = 1*x^6 + 4*x^5 + 6*x^4 + 2*x^3 + 1*x^2 + 2*x + 1 is {N}")

if __name__ == "__main__":
    find_conductor()
