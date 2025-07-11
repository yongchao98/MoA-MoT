import cypari2
import math

def solve_elliptic_curve_problem():
    """
    Solves the problem of finding r + a + b for the given elliptic curve and characters.
    """
    try:
        pari = cypari2.Pari()
        # Set a high precision for the calculations
        pari.default("realprecision", 50)

        # 1. Define the elliptic curve E from its Weierstrass equation
        # y^2 + y = x^3 - x^2 - 10x - 20 corresponds to a_coeffs [0, -1, 1, -10, -20]
        E = pari.ellinit([0, -1, 1, -10, -20])

        # 2. Find the rank r of E(Q)
        # By BSD conjecture, the rank is the order of the zero of L(E,s) at s=1.
        r = int(pari.lfunorderzero(E, 1))

        # 3. Define the primitive cubic Dirichlet characters of conductor 7.
        # The character group mod 7 is cyclic of order 6. Let psi be a generator.
        # The two characters of order 3 are chi1 = psi^2 and chi2 = psi^4.
        # They are primitive as the conductor is prime.
        G = pari.znchargroup(7)
        psi = G[1][0]
        chi1 = psi**2
        chi2 = psi**4

        # 4. Calculate the leading coefficient 'a' for L(E, s, chi1) at s=1
        # First, find the order of the zero at s=1
        k1 = int(pari.lfunorderzero(E, chi1, 1))
        # The leading coefficient is L^(k1)(E, 1, chi1) / k1!
        a = pari.lfun(E, chi1, 1, k1) / math.factorial(k1)

        # 5. Calculate the leading coefficient 'b' for L(E, s, chi2) at s=1
        k2 = int(pari.lfunorderzero(E, chi2, 1))
        b = pari.lfun(E, chi2, 1, k2) / math.factorial(k2)

        # 6. Calculate the final sum r + a + b
        total_sum = r + a + b

        # 7. Print the results
        print(f"The rank of E(Q) is r = {r}")
        print(f"The leading coefficient of L(E, s, chi1) at s=1 is a = {a}")
        print(f"The leading coefficient of L(E, s, chi2) at s=1 is b = {b}")
        print("\nThe final equation is r + a + b:")
        # The Pari objects a and b will be printed in full precision
        print(f"{r} + {a} + {b} = {total_sum}")

        # Round the final result to four decimal places
        result_rounded = round(float(total_sum), 4)
        print(f"\nThe value of r + a + b rounded to four decimal places is: {result_rounded}")

    except ImportError:
        print("Please install the cypari2 library using 'pip install cypari2'")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_elliptic_curve_problem()