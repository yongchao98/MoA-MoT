# This script requires the SageMath environment or the 'sagemath' Python library.
# It can be run directly in a SageMath cell or a Python environment where Sage is installed.

from sage.all import polygen, QQ, HyperellipticCurve

def find_conductor():
    """
    Finds the conductor of the elliptic curve associated with the given hyperelliptic curve.
    """
    # Define the polynomial ring over the rational numbers QQ
    R = polygen(QQ, 'x')
    x = R.gen()

    # Coefficients of the polynomial in the curve equation
    c6, c5, c4, c3, c2, c1, c0 = 1, 4, 6, 2, 1, 2, 1

    # Define the polynomial P(x) from the curve equation y^2 = P(x)
    P = c6*x**6 + c5*x**5 + c4*x**4 + c3*x**3 + c2*x**2 + c1*x + c0

    print(f"The curve is defined by the equation y^2 = {c6}*x^6 + {c5}*x^5 + {c4}*x^4 + {c3}*x^3 + {c2}*x^2 + {c1}*x + {c0}")

    # Create the HyperellipticCurve object
    try:
        C = HyperellipticCurve(P)
    except Exception as e:
        print(f"Could not create the hyperelliptic curve: {e}")
        return

    # Get the Jacobian of the curve
    J = C.jacobian()

    # Decompose the Jacobian into its simple factors up to isogeny.
    # The result is a dictionary where keys are the simple abelian varieties
    # and values are their multiplicities.
    decomposition = J.decomposition()

    if not decomposition:
        print("Could not decompose the Jacobian. The conductor of the Jacobian itself is a larger number.")
        print(f"Jacobian Conductor: {J.conductor()}")
        return

    # The decomposition shows that the Jacobian is isogenous to a product of
    # two copies of a single elliptic curve. We are interested in the conductor
    # of this base elliptic curve.
    
    # Extract the first (and only) simple factor from the decomposition.
    elliptic_curve_factor = list(decomposition.keys())[0]

    # Get the conductor of this elliptic curve factor.
    conductor = elliptic_curve_factor.conductor()

    print(f"The conductor of the associated elliptic curve is: {conductor}")

if __name__ == "__main__":
    find_conductor()