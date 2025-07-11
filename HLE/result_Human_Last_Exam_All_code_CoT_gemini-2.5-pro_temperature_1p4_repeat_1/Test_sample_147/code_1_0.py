import fractions

def solve_optimal_combination():
    """
    Calculates the constants (C, n, m) for the error term of the optimal
    linear combination of Simpson's 1/3 and 3/8 rules.
    """
    # Let the interval [a,b] be divided into 6 subintervals of width k.
    # So, b-a = 6k.

    # Step 1: Error expansion coefficients for the composite rules.
    # The error E = I - Q is of the form A*k^5*f^(4) + B*k^7*f^(6) + ...

    # For Composite Simpson's 1/3 rule (3 applications, step size h=k):
    # From the Euler-Maclaurin formula, the error is:
    # E_1/3 = - (b-a)/180 * h^4 * f^(4) - (b-a)/1512 * h^6 * f^(6)
    # Substituting b-a=6k and h=k:
    # E_1/3 = - (6k)/180 * k^4 * f^(4) - (6k)/1512 * k^6 * f^(6)
    A1 = fractions.Fraction(-6, 180)  # equal to -1/30
    B1 = fractions.Fraction(-6, 1512) # equal to -1/252

    # For Composite Simpson's 3/8 rule (2 applications, step size h=k):
    # The Euler-Maclaurin formula derivation for this rule yields:
    # E_3/8 = -3/40 * k^5 * f^(4) - 53/560 * k^7 * f^(6)
    A2 = fractions.Fraction(-3, 40)
    B2 = fractions.Fraction(-53, 560)

    # Step 2: Find the optimal combination coefficients alpha and beta.
    # Let Q = alpha*Q_1/3 + beta*Q_3/8, with beta = 1 - alpha.
    # The error is E_Q = alpha*E_1/3 + (1-alpha)*E_3/8.
    # We choose alpha to cancel the k^5 * f^(4) term:
    # alpha*A1 + (1-alpha)*A2 = 0  => alpha*(A1 - A2) = -A2
    alpha = A2 / (A2 - A1)
    beta = 1 - alpha

    # Step 3: Calculate the coefficient of the resulting error term.
    # The new error is dominated by the k^7 * f^(6) term:
    # E_Q = (alpha*B1 + beta*B2) * k^7 * f^(6)
    # So, m = 6.
    m = 6
    D_new = alpha * B1 + beta * B2

    # Step 4: Convert to the final error form C*(b-a)^n * f^(m).
    # We have E_Q = D_new * k^7 * f^(6).
    # Substitute k = (b-a)/6:
    # E_Q = D_new * ((b-a)/6)^7 * f^(6)
    # E_Q = (D_new / 6^7) * (b-a)^7 * f^(6)
    # So, n = 7 and C = D_new / 6^7.
    n = 7
    C = D_new / (6**7)

    print(f"The order of the derivative in the error term is m = {m}.")
    print(f"The power of (b-a) in the error term is n = {n}.")
    print(f"The constant C in the error term is the fraction: {C.numerator}/{C.denominator}.")
    
    # Final formatted tuple for the answer
    # This part won't be executed, but shows the final values.
    # return (f"{C.numerator}/{C.denominator}", n, m)

solve_optimal_combination()
<<<('1/4082400', 7, 6)>>>