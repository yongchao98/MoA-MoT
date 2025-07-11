import math

def solve_integral_approximation():
    """
    This script develops and presents an analytical approximation for the integral
    I(eps) = integral from 0 to 15 of 1/(eps + 9x^5 + 5x^6 + 9x^8) dx
    for the small epsilon regime.
    """

    # --- Constants from the problem's integrand ---
    # The term f(x) in the denominator is 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0
    # For small x, the leading term is a*x^p
    p0 = 5.0
    a0 = 9.0

    # --- Derivation of the exponent p ---
    # The exponent in the approximation I(eps) ~ eps^(-p) is given by
    # p = (p0 - 1) / p0
    p_num = p0 - 1.0
    p_den = p0
    p = p_num / p_den

    # --- Derivation of the constant C ---
    # The coefficient C is given by the formula: C = pi / (p0 * a0^(1/p0) * sin(pi/p0))
    try:
        c_val = math.pi / (p0 * (a0**(1.0/p0)) * math.sin(math.pi/p0))
    except (ValueError, ZeroDivisionError) as e:
        print(f"Error calculating C: {e}")
        return

    # --- Presenting the final formula ---
    print("The analytical formula that approximates the integral I(epsilon) for small epsilon is developed as follows:")
    print("1. For small epsilon, the integral's value is dominated by the behavior of the integrand near x=0.")
    print(f"2. The denominator is approximated by its leading term for small x: epsilon + {a0}*x^{p0}.")
    print("3. The integral from 0 to 15.0 is approximated by an integral from 0 to infinity, which can be solved analytically.")
    print("\nThe resulting approximation has the form: I(epsilon) ~= C * epsilon^(-p)\n")

    print("--- Formula Components ---")
    print(f"The exponent is p = ({int(p_num)}.{0}) / ({int(p_den)}.{0})")
    print(f"p = {p}")
    print("")

    print("The constant C is given by the expression:")
    print(f"C = pi / ({p0} * {a0}^(1.0/{p0}) * sin(pi/{p0}))")
    print(f"Numerically, C ~= {c_val:.6f}")
    print("")

    print("--- Final Analytical Formula ---")
    print("The final approximate formula is:")
    print(f"I(epsilon) ~= {c_val:.6f} * epsilon^(-{p})")
    print("\nOr, in a more symbolic form showing the numbers from the problem:")
    print(f"I(epsilon) ~= (pi / ({p0} * {a0}^(1.0/{p0}) * sin(pi/{p0}))) * epsilon^(-{int(p_num)}.{0}/{int(p_den)}.{0})")

solve_integral_approximation()