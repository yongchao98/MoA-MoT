def find_potts_fkg_counterexample():
    """
    This function demonstrates that the local FKG condition for the Potts model fails for q>=3.
    It checks the inequality:
    I{max(a,c)=max(b,d)} + I{min(a,c)=min(b,d)} >= I{a=b} + I{c=d}
    for a, b, c, d in {1, 2, 3}.
    A counterexample corresponds to configurations xi = (a, b) and eta = (c, d)
    on a single edge <u,v>.
    """
    q = 3
    states = range(1, q + 1)
    
    # We will check the known counterexample from the explanation
    a, b, c, d = 2, 2, 1, 3
    
    print("Checking the FKG inequality for the Potts model on a single edge.")
    print(f"Let q = {q}. We will test a specific configuration that provides a counterexample.")
    
    lhs_val_min = min(a, c)
    rhs_val_min = min(b, d)
    lhs_val_max = max(a, c)
    rhs_val_max = max(b, d)
    
    indicator_min = 1 if lhs_val_min == rhs_val_min else 0
    indicator_max = 1 if lhs_val_max == rhs_val_max else 0
    lhs = indicator_min + indicator_max
    
    indicator_a_b = 1 if a == b else 0
    indicator_c_d = 1 if c == d else 0
    rhs = indicator_a_b + indicator_c_d

    print("\nLet xi(u) = a, xi(v) = b for one configuration xi.")
    print("Let eta(u) = c, eta(v) = d for another configuration eta.")
    print(f"Consider the values: a = {a}, b = {b}, c = {c}, d = {d}.")

    print("\nThe inequality to check is:")
    print("I{max(a,c) = max(b,d)} + I{min(a,c) = min(b,d)} >= I{a=b} + I{c=d}\n")
    
    # Detailed LHS calculation
    print("Calculating the Left-Hand Side (LHS):")
    print(f"max(a,c) = max({a},{c}) = {lhs_val_max}")
    print(f"max(b,d) = max({b},{d}) = {rhs_val_max}")
    print(f"Is max(a,c) = max(b,d)? {lhs_val_max} = {rhs_val_max}? {'Yes' if indicator_max == 1 else 'No'}. Indicator = {indicator_max}.")
    print(f"min(a,c) = min({a},{c}) = {lhs_val_min}")
    print(f"min(b,d) = min({b},{d}) = {rhs_val_min}")
    print(f"Is min(a,c) = min(b,d)? {lhs_val_min} = {rhs_val_min}? {'Yes' if indicator_min == 1 else 'No'}. Indicator = {indicator_min}.")
    print(f"LHS = {indicator_max} + {indicator_min} = {lhs}")
    
    # Detailed RHS calculation
    print("\nCalculating the Right-Hand Side (RHS):")
    print(f"Is a = b? {a} = {b}? {'Yes' if indicator_a_b == 1 else 'No'}. Indicator = {indicator_a_b}.")
    print(f"Is c = d? {c} = {d}? {'Yes' if indicator_c_d == 1 else 'No'}. Indicator = {indicator_c_d}.")
    print(f"RHS = {indicator_a_b} + {indicator_c_d} = {rhs}")
    
    # Final check
    print(f"\nComparing LHS and RHS:")
    print(f"Is {lhs} >= {rhs}? This is {'True' if lhs >= rhs else 'False'}.")
    print("\nThe inequality fails, providing a counterexample.")

find_potts_fkg_counterexample()