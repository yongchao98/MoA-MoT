import math

def compute_generalized_markov_cf(p, q):
    """
    Computes the continued fraction for the generalized Markov number m_{p/q}.
    """
    # Step 1 & 2: Calculate the number m_{p/q} = sqrt(Delta)/q
    delta_orig = 9 * p**2 - 8 * p * q + 4 * q**2
    
    # Step 3: Prepare for the continued fraction algorithm
    # We use D = delta_orig * q^2 to ensure integer arithmetic
    D = delta_orig * q**2
    sqrt_D = math.isqrt(D)
    
    # Initial state for the algorithm: x = (alpha_0 + sqrt(D)) / beta_0
    alpha = 0
    beta = q**2
    
    # Check if the initial state satisfies the algorithm's divisibility condition
    if (D - alpha**2) % beta != 0:
        print(f"Error: The initial state for p={p}, q={q} does not satisfy the algorithm's requirements.")
        return

    # Step 4 & 5: Implement the algorithm and detect periodicity
    coeffs = []
    seen_states = {}
    
    while (alpha, beta) not in seen_states:
        seen_states[(alpha, beta)] = len(coeffs)
        
        # Calculate the next coefficient a_k
        a = (alpha + sqrt_D) // beta
        coeffs.append(a)
        
        # Calculate the next state (alpha_{k+1}, beta_{k+1})
        alpha_next = a * beta - alpha
        beta_next = (D - alpha_next**2) // beta
        
        alpha = alpha_next
        beta = beta_next

    # The loop terminates when a state is repeated. Find where the repetition starts.
    start_of_period = seen_states[(alpha, beta)]
    
    integer_part = coeffs[:start_of_period]
    repeating_part = coeffs[start_of_period:]

    # Step 6: Output the results
    print(f"The continued fraction for the generalized Markov number m_{p}/{q} is computed.")
    
    # Combine integer and repeating parts for the standard notation string
    integer_str = ", ".join(map(str, integer_part[1:]))
    repeating_str = ", ".join(map(str, repeating_part))
    
    if integer_part:
        if len(integer_part) > 1:
            cf_notation = f"[{integer_part[0]}; {integer_str}, ({repeating_str})]"
        else:
            cf_notation = f"[{integer_part[0]}; ({repeating_str})]"
    else: # Should not happen for this problem
        cf_notation = f"[0; ({repeating_str})]"
    
    print(f"\nThe standard notation is: {cf_notation}\n")

    print("This means the numbers in the final equation are:")
    if integer_part:
        print(f"The integer part is a_0 = {integer_part[0]}")
    
    print("The repeating block of coefficients is:")
    print(', '.join(map(str, repeating_part)))


# Compute for m_{4/7}
p, q = 4, 7
compute_generalized_markov_cf(p, q)