def check_potts_fkg_condition():
    """
    Checks the Holley's criterion for the Potts model on a single edge.
    The positive correlation property holds if and only if the following
    inequality is true for all spin assignments a,b,c,d.
    Inequality: I(max(a,c)==max(b,d)) + I(min(a,c)==min(b,d)) >= I(a==b) + I(c==d)
    
    This script searches for a counterexample for q=3.
    """
    q = 3
    print(f"Verifying the positive correlation condition for the Potts model.")
    print(f"The condition boils down to checking an inequality for all spin values on a single edge.")
    print(f"We search for a counterexample with q = {q} states.")
    
    # Iterate through all combinations of spin values a, b, c, d from {1, ..., q}
    for a in range(1, q + 1):
        for b in range(1, q + 1):
            for c in range(1, q + 1):
                for d in range(1, q + 1):
                    # Using the fact that bools evaluate to 1 (True) or 0 (False)
                    lhs = (max(a, c) == max(b, d)) + (min(a, c) == min(b, d))
                    rhs = (a == b) + (c == d)
                    
                    if lhs < rhs:
                        print("\nA counterexample to the inequality has been found:")
                        print(f"Let a={a}, b={b}, c={c}, d={d}.")

                        # Store values for detailed print
                        val_max_ac = max(a,c)
                        val_max_bd = max(b,d)
                        val_min_ac = min(a,c)
                        val_min_bd = min(b,d)
                        lhs_term1 = 1 if (val_max_ac == val_max_bd) else 0
                        lhs_term2 = 1 if (val_min_ac == val_min_bd) else 0
                        rhs_term1 = 1 if (a == b) else 0
                        rhs_term2 = 1 if (c == d) else 0
                        
                        print("\nEvaluating the inequality:")
                        print(f"I(max({a},{c}) == max({b},{d})) + I(min({a},{c}) == min({b},{d})) >= I({a}=={b}) + I({c}=={d})")
                        print(f"I({val_max_ac} == {val_max_bd}) + I({val_min_ac} == {val_min_bd}) >= {rhs_term1} + {rhs_term2}")
                        print(f"{lhs_term1} + {lhs_term2} >= {rhs}")
                        print(f"{lhs} >= {rhs}")
                        print(f"This is FALSE.")
                        
                        print("\nThis single failure implies the positive correlation property does not hold for q>=3 on any graph with one or more edges.")
                        print("A graph has at least one edge if its maximum degree is d >= 1.")
                        print("The property does hold for d=0 (a graph with a single vertex).")
                        print("Therefore, the largest integer d for which the property holds is 0.")
                        return

# Run the function to find the counterexample and print the conclusion.
check_potts_fkg_condition()
