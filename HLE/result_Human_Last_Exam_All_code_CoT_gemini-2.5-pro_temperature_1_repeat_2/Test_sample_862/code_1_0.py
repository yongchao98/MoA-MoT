import math
from fractions import Fraction

def solve_for_C():
    """
    This function calculates the smallest possible constant C by solving the
    underlying eigenvalue problem.
    """
    
    # The problem is to find the smallest constant C such that
    # integral(a * f^2) <= C * integral(a * (f')^2)
    # for a(x) in [1, 3] and integral(a * f) = 0.
    # This constant C is the supremum of the Rayleigh quotient, which is
    # C = sup_a(x) (1 / lambda_1(a)), where lambda_1(a) is the first non-zero
    # eigenvalue of the Sturm-Liouville operator L(f) = -(1/a)(a*f')'
    # with periodic boundary conditions.
    # This is equivalent to C = 1 / inf_a(x) (lambda_1(a)).
    # The infimum of lambda_1 is achieved for a bang-bang control function a(x)
    # which takes values m=1 and M=3. The associated eigenfunction can be
    # either even or odd. We must find the lowest eigenvalue from both cases.

    # Let lambda = k^2. The analysis leads to solving a transcendental equation
    # for k, which is minimized with respect to the switching point of a(x).
    # This leads to two systems of equations, one for the even case and one for the odd case.

    # Define the bounds for a(x)
    m = 1
    M = 3

    print("Solving for the minimal eigenvalue lambda_1 = k^2.")
    print(f"The weight function a(x) is in the range [{m}, {M}].\n")

    # --- Even eigenfunction case ---
    # The analysis leads to the following system:
    # sin^2(u) = m / (M + m)
    # sin^2(v) = M / (M + m)
    # where u = k*l, v = k*(pi-l), u in (0, pi/2), v in (pi/2, pi).
    # The parameter l is the switching point for the bang-bang function a(x).
    print("--- Case 1: Even eigenfunction ---")

    # From sin^2(u) = 1/4 and u in (0, pi/2), we have u = pi/6.
    # From sin^2(v) = 3/4 and v in (pi/2, pi), we have v = 2*pi/3.
    u_even = Fraction(1, 6) # in units of pi
    v_even = Fraction(2, 3) # in units of pi

    # The eigenvalue k_even is given by k_even * pi = u + v
    k_even = u_even + v_even
    lambda_even = k_even**2

    print(f"Found solution: u = ({u_even})*pi, v = ({v_even})*pi")
    print(f"This gives k_even = (u+v)/pi = {k_even}")
    print(f"The corresponding eigenvalue is lambda_even = k_even^2 = ({k_even.numerator}/{k_even.denominator})^2 = {lambda_even.numerator}/{lambda_even.denominator}\n")


    # --- Odd eigenfunction case ---
    # The analysis leads to a similar system:
    # sin^2(u) = M / (M + m)
    # sin^2(v) = m / (M + m)
    # where u = k*l, v = k*(pi-l), u in (0, pi/2), v in (pi/2, pi).
    print("--- Case 2: Odd eigenfunction ---")

    # From sin^2(u) = 3/4 and u in (0, pi/2), we have u = pi/3.
    # From sin^2(v) = 1/4 and v in (pi/2, pi), we have v = 5*pi/6.
    u_odd = Fraction(1, 3) # in units of pi
    v_odd = Fraction(5, 6) # in units of pi

    # The eigenvalue k_odd is given by k_odd * pi = u + v
    k_odd = u_odd + v_odd
    lambda_odd = k_odd**2

    print(f"Found solution: u = ({u_odd})*pi, v = ({v_odd})*pi")
    print(f"This gives k_odd = (u+v)/pi = {k_odd}")
    print(f"The corresponding eigenvalue is lambda_odd = k_odd^2 = ({k_odd.numerator}/{k_odd.denominator})^2 = {lambda_odd.numerator}/{lambda_odd.denominator}\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    min_lambda = min(lambda_even, lambda_odd)
    print(f"The minimum first non-zero eigenvalue is min({lambda_even.numerator}/{lambda_even.denominator}, {lambda_odd.numerator}/{lambda_odd.denominator}) = {min_lambda.numerator}/{min_lambda.denominator}.")

    C = 1 / min_lambda
    print(f"The smallest possible constant C is 1 / lambda_min = 1 / ({min_lambda.numerator}/{min_lambda.denominator}) = {C.numerator}/{C.denominator}.")

    # Output the final answer
    print("\nThe final equation is C = num/den. The numbers are:")
    print(f"Numerator: {C.numerator}")
    print(f"Denominator: {C.denominator}")

solve_for_C()