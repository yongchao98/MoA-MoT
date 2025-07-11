import math

def compute_markov_cf():
    """
    Computes and prints the continued fraction associated with the generalized
    Markov number m_{4/7}.
    """
    p = 4
    q = 7

    # 1. Define the quadratic irrational alpha = (C - sqrt(D)) / B
    D = 9 * q**2 - 4 * p**2
    C = 3 * q
    B = 2 * p

    # 2. Implement the continued fraction algorithm.
    # The terms xi_i for i >= 1 will be of the form (m + sqrt(D)) / d.
    
    sqrt_D_int = math.isqrt(D)
    
    # First term, a_0
    # alpha = (21 - sqrt(377)) / 8 approx (21 - 19.4) / 8 = 0.2
    a0 = (C - sqrt_D_int) // B
    
    # Setup for the periodic part. The first term of the periodic part is
    # xi_1 = 1 / (alpha - a0) = B / (C - sqrt(D)) = (C + sqrt(D)) / ((C^2-D)/B)
    m = C
    d = (C**2 - D) // B
    
    # Store coefficients and states to find the cycle
    all_coeffs = [a0]
    seen_states = {}
    
    while (m, d) not in seen_states:
        seen_states[(m, d)] = len(all_coeffs)
        
        # Calculate next coefficient a_i = floor((m + sqrt(D)) / d)
        a = (m + sqrt_D_int) // d
        all_coeffs.append(a)
        
        # Calculate next state (m_{i+1}, d_{i+1})
        m_next = a * d - m
        d_next = (D - m_next**2) // d
        
        m = m_next
        d = d_next

    # 3. Format the output string.
    cycle_start_index = seen_states[(m, d)]
    
    non_periodic_part = all_coeffs[:cycle_start_index]
    repeating_part = all_coeffs[cycle_start_index:]
    
    head = non_periodic_part[0]
    # The rest of the non-periodic part (if any)
    tail_str = ", ".join(map(str, non_periodic_part[1:]))
    repeating_str = ", ".join(map(str, repeating_part))

    # Construct the final continued fraction string
    cf_str = f"[{head};"
    if tail_str:
        cf_str += f" {tail_str},"
    cf_str += f" ({repeating_str})]"

    # 4. Print the final equation with all its numbers.
    print(f"The final equation for the continued fraction associated with m_{p}/{q} is:")
    print(f"({C} - sqrt({D})) / {B} = {cf_str}")

compute_markov_cf()
<<<[0; 5, (19, 4, 1, 4)]>>>