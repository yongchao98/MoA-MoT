import math

def solve_task():
    """
    Computes the continued fraction for the generalized Markov number m_{4/7}.
    """
    p, q = 4, 7

    # 1. Define the quadratic equation q*x^2 - (3q - 2p)*x + p = 0
    # Coefficients of the quadratic equation ax^2 + bx + c = 0
    a = q
    b = -(3 * q - 2 * p)
    c = p

    # 2. Solve for the larger root, which is m_{p/q}
    # x = (-b +/- sqrt(b^2 - 4ac)) / 2a
    discriminant = b**2 - 4 * a * c
    # Numerator of the root will be of the form P0 + sqrt(D)
    # The number is (P0_numerator / P0_denominator) + (sqrt(D) / Q0)
    # We want it in the form (P0 + sqrt(D)) / Q0
    # The larger root is (-b + sqrt(discriminant)) / (2a)
    P0 = -b
    Q0 = 2 * a
    D = discriminant
    
    # Simplify the fraction by dividing by GCD if possible
    common_divisor = math.gcd(P0, Q0)
    # We can't simplify the sqrt(D) term this way unless D is a perfect square
    # and we know it's not. The algorithm for continued fractions handles this.

    # 3. Compute the continued fraction using integer arithmetic
    # The algorithm works for numbers of the form (P + sqrt(D)) / Q
    
    coeffs = []
    seen_states = {}
    
    # Initial state
    Pi = P0
    Qi = Q0

    sqrt_D_floor = math.isqrt(D)

    while (Pi, Qi) not in seen_states:
        seen_states[(Pi, Qi)] = len(coeffs)
        
        # Calculate the integer part
        ai = (Pi + sqrt_D_floor) // Qi
        coeffs.append(ai)
        
        # Calculate the next state (Pi+1, Qi+1)
        # x_i+1 = 1 / (x_i - a_i)
        # Let x_i = (Pi + sqrt(D)) / Qi
        # x_i - a_i = (Pi - a_i*Qi + sqrt(D)) / Qi
        # 1 / (x_i - a_i) = Qi / (Pi - a_i*Qi + sqrt(D))
        # Rationalize: Qi * (sqrt(D) - (Pi - a_i*Qi)) / (D - (Pi - a_i*Qi)^2)
        # This gives Pi+1 = a_i*Qi - Pi
        # And Qi+1 = (D - Pi+1^2) / Qi
        
        Pi_next = ai * Qi - Pi
        Qi_next = (D - Pi_next**2) // Qi
        
        Pi, Qi = Pi_next, Qi_next

    # 4. Format the output string
    start_of_period = seen_states[(Pi, Qi)]
    
    integer_part = coeffs[0]
    non_periodic_part = coeffs[1:start_of_period]
    periodic_part = coeffs[start_of_period:]

    # Build the continued fraction string representation
    cf_parts = []
    if non_periodic_part:
        cf_parts.append(", ".join(map(str, non_periodic_part)))
    if periodic_part:
        cf_parts.append(f"({', '.join(map(str, periodic_part))})")
    
    cf_string = f"[{integer_part}; {', '.join(cf_parts)}]"

    # Final result string
    result = (f"The generalized Markov number m_{p}/{q} is the larger root of "
              f"{a}x^2 + ({b})x + {c} = 0.\n"
              f"m_{p}/{q} = ({P0} + \u221A{D}) / {Q0}\n"
              f"The associated continued fraction is: {cf_string}")

    print(result)

solve_task()
<<<[1; 2, (7, 3, 1, 1, 1, 3)]>>>