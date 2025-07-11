import math

def solve_generalized_markov_cf(p, q):
    """
    Computes the continued fraction associated with the generalized Markov number m_{p/q}.

    This corresponds to the larger root of the quadratic equation:
    q*x^2 - (3*q - 2*p)*x + p = 0
    """
    # 1. Define coefficients of the quadratic equation ax^2 + bx + c = 0
    a = q
    b = -(3 * q - 2 * p)
    c = p

    print(f"For p/q = {p}/{q}, the associated quadratic equation is:")
    print(f"{a}x^2 + ({b})x + ({c}) = 0")
    print("-" * 30)

    # 2. Find the components of the larger root (P + sqrt(D)) / Q
    # The root is (-b + sqrt(b^2 - 4ac)) / (2a)
    D = b**2 - 4 * a * c
    P0 = -b
    Q0 = 2 * a

    # Simplify the fraction by finding the greatest common divisor
    common_divisor = math.gcd(math.gcd(P0, Q0), int(math.sqrt(D)) if D > 0 and math.isqrt(D)**2 == D else 1)
    # In this problem, D=57, which is not a perfect square, so common_divisor is just gcd(P0, Q0)
    common_divisor = math.gcd(P0, Q0)
    #P0 //= common_divisor
    #Q0 //= common_divisor
    #D //= (common_divisor**2) # This is only for perfect square D

    print(f"The larger root is ({P0} + sqrt({D})) / {Q0}")
    print("-" * 30)

    # 3. Compute the continued fraction using integer arithmetic
    coeffs = []
    seen_states = {}
    
    Pi = P0
    Qi = Q0
    sqrt_D_int = math.isqrt(D)

    i = 0
    while True:
        state = (Pi, Qi)
        if state in seen_states:
            start_of_period = seen_states[state]
            break
        
        seen_states[state] = i
        
        # a_i = floor((Pi + sqrt(D)) / Qi)
        ai = (Pi + sqrt_D_int) // Qi
        coeffs.append(ai)
        
        # Update P and Q for the next iteration
        # P_{i+1} = a_i * Q_i - P_i
        # Q_{i+1} = (D - P_{i+1}^2) / Q_i
        P_next = ai * Qi - Pi
        Q_next = (D - P_next**2) // Qi
        
        Pi = P_next
        Qi = Q_next
        i += 1

    # 4. Format and print the result
    non_periodic_part = coeffs[:start_of_period]
    periodic_part = coeffs[start_of_period:]

    print("The continued fraction is:")
    
    # Print the equation part by part
    non_periodic_str = ', '.join(map(str, non_periodic_part[1:]))
    periodic_str = ', '.join(map(str, periodic_part))
    
    result_str = f"m_{p}/{q} = [{non_periodic_part[0]}; "
    if non_periodic_part[1:]:
        result_str += f"{non_periodic_str}, "
    result_str += f"({periodic_str})]"

    print(result_str)
    
    print("\nWhere the numbers in the final equation are:")
    print(f"Integer part: {non_periodic_part[0]}")
    if len(non_periodic_part) > 1:
        print(f"Non-repeating part (after integer part): {', '.join(map(str, non_periodic_part[1:]))}")
    print(f"Repeating part: {', '.join(map(str, periodic_part))}")
    
    return f"[{non_periodic_part[0]}; {non_periodic_str}, ({periodic_str})]"


if __name__ == '__main__':
    p_val = 4
    q_val = 7
    solve_generalized_markov_cf(p_val, q_val)
