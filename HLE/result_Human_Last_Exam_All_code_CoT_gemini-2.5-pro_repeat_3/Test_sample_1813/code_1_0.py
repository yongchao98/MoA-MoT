import math

def compute_continued_fraction_for_markov():
    """
    Computes the continued fraction for the generalized Markov number m_{4/7}.
    """
    r, s = 4, 7

    # Calculate the components of the number alpha = (p0 + sqrt(d0)) / q0
    p0 = 3 * s - 2 * r
    d0 = 9 * s**2 - 4
    q0 = 2 * s

    print(f"The number is alpha = ({p0} + sqrt({d0})) / {q0}")

    # Step 1: Compute a0
    sqrt_d0 = math.sqrt(d0)
    a0 = math.floor((p0 + sqrt_d0) / q0)
    
    coeffs = [a0]
    
    # Step 2: Transform the remainder to a form suitable for the PQa algorithm.
    # The remainder is xi_1 = 1 / (alpha - a0).
    # After rationalizing and rearranging, we get xi_1 = (P1 + sqrt(D1)) / Q1
    # P1 = a0*q0^2 - p0*q0
    # D1 = d0*q0^2
    # Q1 = d0 - (p0 - a0*q0)^2
    # To avoid large numbers, we simplify this transformation.
    # xi_1 = q0 / (p0 - a0*q0 + sqrt(d0))
    # Let num = p0 - a0*q0 = 13 - 2*14 = -15
    # xi_1 = q0 * (-num + sqrt(d0)) / (d0 - num^2)
    # xi_1 = (14 * 15 + 14 * sqrt(437)) / (437 - 225)
    # xi_1 = (210 + 14 * sqrt(437)) / 212
    # xi_1 = (105 + 7 * sqrt(437)) / 106
    # xi_1 = (105 + sqrt(49*437)) / 106
    # xi_1 = (105 + sqrt(21413)) / 106
    P = 105
    D = 21413
    Q = 106

    sqrt_D = math.sqrt(D)
    
    # Keep track of (P, Q) pairs to detect the cycle
    seen_states = {}
    
    # Step 3: Run the PQa algorithm
    while (P, Q) not in seen_states:
        seen_states[(P, Q)] = len(coeffs)
        
        a = math.floor((P + sqrt_D) / Q)
        coeffs.append(a)
        
        P_next = a * Q - P
        Q_next = (D - P_next**2) // Q
        
        P = P_next
        Q = Q_next

    # The cycle starts at the index stored in seen_states[(P,Q)]
    start_of_period = seen_states[(P, Q)]
    non_repeating_part = coeffs[:start_of_period]
    repeating_part = coeffs[start_of_period:]

    a0_str = str(non_repeating_part[0])
    repeating_str = ", ".join(map(str, repeating_part))

    print("\nThe continued fraction is:")
    # Using overline for repeating part is standard but hard in plain text.
    # We will show the equation with the coefficients.
    print(f"({p0} + sqrt({d0})) / {q0} = [{a0_str}; {repeating_str}, ...]")
    print("\nEach number in the final equation:")
    print(f"Integer part a_0: {a0_str}")
    print(f"Repeating part: {repeating_part}")

compute_continued_fraction_for_markov()
<<<[2; 2, 2, 1, 2, 3, 1, 1, 2]>>>