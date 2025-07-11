def N(coeffs):
    """
    Computes the numerator of a continued fraction defined by a list of coefficients.
    N[a_1, ..., a_k] is calculated using the recurrence p_n = a_n * p_{n-1} + p_{n-2},
    with initial values p_0 = 1, p_{-1} = 0.
    """
    if not isinstance(coeffs, list):
        raise TypeError("Input must be a list of coefficients.")
    
    p_prev2 = 1  # Corresponds to p_0
    p_prev1 = 0  # Dummy value to start, will be updated to p_1 first
    
    if len(coeffs) == 0:
        return p_prev2 # N([]) = p_0 = 1
        
    p_prev1 = coeffs[0] # p_1 = a_1*p_0 + p_{-1} if p_0=1,p-1=0? No. Let's use simpler logic.
    
    # Let's re-verify the base cases for the loop logic
    # p_n = a_n p_{n-1} + p_{n-2}
    # Base cases: p_0=1, p_1=a_1
    
    if len(coeffs) == 0:
      return 1 # p_0
    
    p_prev_val = 1 # p_0
    p_curr_val = coeffs[0] # p_1
    
    if len(coeffs) == 1:
      return p_curr_val
      
    for i in range(1, len(coeffs)):
        p_next_val = coeffs[i] * p_curr_val + p_prev_val
        p_prev_val = p_curr_val
        p_curr_val = p_next_val
        
    return p_curr_val

def solve_for_ck():
    """
    Solves for c_k for a sample case and verifies the derived formula.
    """
    # Example case: k=4, a_i = [1, 2, 3, 2]
    k = 4
    a = [1, 2, 3, 2] # a_1, a_2, a_3, a_4
    
    print(f"Solving for k={k} with coefficients a = {a}")

    # Construct sequences for the main equation
    # LHS sequence: [a_2, ..., a_{k-1}, a_k+1, a_k, ..., a_1]
    seq_LHS = a[1:k-1] + [a[k-1] + 1] + a[::-1]
    
    # RHS sequence: [a_1, ..., a_k, a_k, ..., a_2]
    seq_RHS = a[:k] + a[1:k][::-1]

    # Calculate LHS and RHS numerators
    val_LHS = N(seq_LHS)
    val_RHS = N(seq_RHS)

    # Calculate c_k from the original equation
    c_k_from_eq = val_LHS - val_RHS

    # Calculate c_k using our derived formula
    # c_k = N[a_1, ..., a_k] * N[a_2, ..., a_{k-1}]
    term1 = N(a[:k])
    term2 = N(a[1:k-1])
    c_k_from_formula = term1 * term2

    # Print the results
    print("\nVerifying the equation:")
    print(f"N{seq_LHS} = {val_LHS}")
    print(f"N{seq_RHS} = {val_RHS}")
    print(f"From the equation, c_k = {val_LHS} - {val_RHS} = {c_k_from_eq}")
    print(f"The equation with numbers is: {val_LHS} = {val_RHS} + {c_k_from_eq}\n")
    
    print("Verifying the derived formula for c_k:")
    print(f"c_k = N[a_1,...,a_k] * N[a_2,...,a_k-1]")
    print(f"N{a[:k]} = {term1}")
    print(f"N{a[1:k-1]} = {term2}")
    print(f"From the formula, c_k = {term1} * {term2} = {c_k_from_formula}")

    if c_k_from_eq == c_k_from_formula:
        print("\nThe results match. The formula for c_k is correct.")
    else:
        print("\nThere is a discrepancy in the results.")

solve_for_ck()
