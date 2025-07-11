import math

def compute_markov_continued_fraction(p, q):
    """
    Computes the continued fraction associated with the generalized Markov number m_{p/q}.

    This is the continued fraction of the larger root of the equation:
    q*x^2 - (3*q - p)*x - q = 0
    """
    # 1. Set up the quadratic equation ax^2 + bx + c = 0
    a = q
    b = -(3 * q - p)
    c = -q

    print(f"The associated quadratic equation for m_{p}/{q} is:")
    # The instruction says to output each number in the final equation.
    # Here are the coefficients of the quadratic equation:
    print(f"{a}x^2 + ({b})x + ({c}) = 0")
    print(f"Which simplifies to: {a}x^2 - {-b}x - {-c} = 0")
    print("-" * 30)

    # 2. Find the larger root using the quadratic formula
    # x = (-b +/- sqrt(b^2 - 4ac)) / 2a
    discriminant = b**2 - 4 * a * c
    numerator_const = -b
    denominator = 2 * a
    
    # We work with the form (P + sqrt(D)) / Q
    # Initial P, Q, D for the algorithm
    P_0 = numerator_const
    Q_0 = denominator
    D = discriminant

    print(f"The larger root is x = ({P_0} + sqrt({D})) / {Q_0}")
    print("-" * 30)
    
    # 3. Compute the continued fraction
    terms = []
    seen_states = {}
    
    P = P_0
    Q = Q_0
    isqrt_D = math.isqrt(D)

    while True:
        state = (P, Q)
        if state in seen_states:
            start_of_cycle = seen_states[state]
            non_repeating_part = terms[:start_of_cycle]
            repeating_part = terms[start_of_cycle:]
            break
        
        seen_states[state] = len(terms)
        
        # Calculate the integer part
        term = (P + isqrt_D) // Q
        terms.append(term)
        
        # Update P and Q for the next iteration
        # P_next = a*Q - P
        # Q_next = (D - P_next^2) / Q
        P_next = term * Q - P
        Q_next = (D - P_next**2) // Q
        
        P = P_next
        Q = Q_next

    print("The continued fraction expansion is:")
    # The instruction says to output each number in the final equation.
    # Interpreting this as the coefficients of the continued fraction:
    if not non_repeating_part:
        # Purely periodic
        repeating_str = ", ".join(map(str, repeating_part))
        print(f"[{repeating_str}] (periodic)")
        print(f"The numbers in the repeating block are: {repeating_part}")
    else:
        non_repeating_str = ", ".join(map(str, non_repeating_part))
        repeating_str = ", ".join(map(str, repeating_part))
        print(f"[{non_repeating_str}; {repeating_str}] (the latter part is periodic)")
        print(f"The number in the non-repeating part is: {non_repeating_part}")
        print(f"The numbers in the repeating block are: {repeating_part}")

# For the generalized Markov number m_{4/7}
p_val = 4
q_val = 7
compute_markov_continued_fraction(p_val, q_val)
