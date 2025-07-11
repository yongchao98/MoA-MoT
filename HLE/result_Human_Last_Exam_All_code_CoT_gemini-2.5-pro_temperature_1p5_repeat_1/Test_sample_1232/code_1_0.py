def verify_potts_fkg_condition():
    """
    This script tests the local condition required for the Potts model
    to have the positive correlations property (FKG inequality).

    The condition is checked for different numbers of states, q.
    The property holds if and only if for every edge, the following inequality
    is true for all spin choices a, b, c, d from {1, ..., q}:
    I(max(a,c)=max(b,d)) + I(min(a,c)=min(b,d)) >= I(a=b) + I(c=d)
    where I() is the indicator function.
    """
    
    print("Verifying the local FKG condition for the Potts model.")
    
    # --- Check for q=2 (Ising model equivalent) ---
    q = 2
    holds_for_q2 = True
    for a in range(1, q + 1):
        for b in range(1, q + 1):
            for c in range(1, q + 1):
                for d in range(1, q + 1):
                    lhs = (1 if max(a, c) == max(b, d) else 0) + \
                          (1 if min(a, c) == min(b, d) else 0)
                    rhs = (1 if a == b else 0) + \
                          (1 if c == d else 0)
                    if lhs < rhs:
                        holds_for_q2 = False
                        break
    
    print(f"\nFor q = {q}, does the condition hold? {holds_for_q2}")
    
    # --- Check for q=3 ---
    q = 3
    found_counterexample_q3 = False
    for a in range(1, q + 1):
        if found_counterexample_q3: break
        for b in range(1, q + 1):
            if found_counterexample_q3: break
            for c in range(1, q + 1):
                if found_counterexample_q3: break
                for d in range(1, q + 1):
                    lhs = (1 if max(a, c) == max(b, d) else 0) + \
                          (1 if min(a, c) == min(b, d) else 0)
                    rhs = (1 if a == b else 0) + \
                          (1 if c == d else 0)
                    if lhs < rhs:
                        print(f"\nFor q = {q}, the condition FAILS.")
                        print("A counterexample is found:")
                        print(f"Let spins be assigned as: a={a}, b={b}, c={c}, d={d}")
                        print("\nThe inequality is: I(max(a,c)=max(b,d)) + I(min(a,c)=min(b,d)) >= I(a=b) + I(c=d)")
                        print("Plugging in the numbers:")
                        print(f"RHS = I({a}={b}) + I({c}={d}) = {1 if a==b else 0} + {1 if c==d else 0} = {rhs}")
                        print(f"LHS = I(max({a},{c})=max({b},{d})) + I(min({a},{c})=min({b},{d}))")
                        print(f"    = I({max(a,c)}={max(b,d)}) + I({min(a,c)}={min(b,d)})")
                        print(f"    = {1 if max(a,c)==max(b,d) else 0} + {1 if min(a,c)==min(b,d) else 0} = {lhs}")
                        print(f"\nThe resulting inequality is {lhs} >= {rhs}, which is FALSE.")
                        found_counterexample_q3 = True
                        break

    print("\nConclusion from the script:")
    print("The property fails for q=3. Since the question requires it to hold for ALL q >= 2,")
    print("it must fail for any graph with at least one edge.")
    print("A graph can only be guaranteed to have no edges if its maximum degree d is 0.")
    print("Therefore, the largest value of d for which the statement holds is 0.")

verify_potts_fkg_condition()