import cypari2
import cmath

def solve_elliptic_curve_problem():
    """
    Solves the problem of finding r + a + b.
    """
    try:
        # Initialize the PARI/GP environment with high precision
        pari = cypari2.Pari()
        pari.default("realprecision", 50)

        # 1. Define the elliptic curve E from its Weierstrass equation
        # y^2 + y = x^3 - x^2 - 10x - 20
        # The coefficients are [a1, a2, a3, a4, a6]
        E = pari.ellinit([0, -1, 1, -10, -20])

        # 2. Compute the rank r of E(Q)
        # The ellrank function returns a vector; the rank is the first component.
        r_pari = pari.ellrank(E)
        r = int(r_pari[0])

        # 3. Find the two primitive cubic Dirichlet characters modulo 7
        G = pari.chargroup(7)
        cubic_chars = [g[0] for g in G if g[1] == 3]
        chi1 = cubic_chars[0]
        chi2 = cubic_chars[1]

        # 4. As explained in the plan, we assume the order of vanishing is 0.
        # The leading coefficients a and b are the L-values at s=1.
        # We use the lfun function to compute the value of the twisted L-series.
        a_pari = pari.lfun(E, chi1, 1)
        b_pari = pari.lfun(E, chi2, 1)

        # Convert the PARI/GP complex numbers to standard Python complex numbers
        a = complex(a_pari)
        b = complex(b_pari)

        # 5. Calculate the final sum r + a + b
        # Since b is the complex conjugate of a, their sum is 2 * Re(a).
        result = r + a + b

        # 6. Output the components of the final equation and the result
        print(f"The rank of the elliptic curve is r = {r}")
        print(f"The leading coefficient a is {a}")
        print(f"The leading coefficient b is {b}")
        print("\nThe final equation is r + a + b:")
        print(f"{r} + ({a}) + ({b}) = {result.real}")
        
        print(f"\nThe final value rounded to four decimal places is: {result.real:.4f}")

    except ImportError:
        print("This script requires the 'cypari2' library.")
        print("Please install it using: pip install cypari2")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure your PARI/GP library is up to date.")

solve_elliptic_curve_problem()