def check_fkg_condition_for_potts(q_max):
    """
    Checks the FKG lattice condition for the Potts model for q from 2 to q_max.

    The condition is checked for a single edge. It must hold for all
    spin assignments a, b, c, d on the endpoints of the edge.
    Let a = xi(x), b = xi(y) and c = eta(x), d = eta(y).
    The inequality to check is:
    I(max(a,c)==max(b,d)) + I(min(a,c)==min(b,d)) >= I(a==b) + I(c==d)
    
    A counterexample to this inequality for a given q proves that the
    positive correlations property does not hold in general for the q-state
    Potts model on graphs with at least one edge.
    """
    for q in range(2, q_max + 1):
        print(f"--- Checking for q = {q} ---")
        spins = range(1, q + 1)
        found_counterexample = False
        
        # Iterate through all possible spin assignments for two configurations on two sites
        for a in spins:
            for b in spins:
                for c in spins:
                    for d in spins:
                        # Left Hand Side of the inequality
                        val_max_ac = max(a, c)
                        val_max_bd = max(b, d)
                        val_min_ac = min(a, c)
                        val_min_bd = min(b, d)
                        lhs = (1 if val_max_ac == val_max_bd else 0) + \
                              (1 if val_min_ac == val_min_bd else 0)
                        
                        # Right Hand Side of the inequality
                        rhs = (1 if a == b else 0) + \
                              (1 if c == d else 0)
                        
                        if lhs < rhs:
                            print("Found a counterexample to the FKG lattice condition:")
                            print(f"  Let the two sites be x and y.")
                            print(f"  Configuration xi: xi(x) = {a}, xi(y) = {b}")
                            print(f"  Configuration eta: eta(x) = {c}, eta(y) = {d}")
                            print("\n  Let's check the inequality for the edge <x,y>:")
                            print(f"  LHS = I(max({a},{c})==max({b},{d})) + I(min({a},{c})==min({b},{d}))")
                            print(f"      = I({val_max_ac}=={val_max_bd}) + I({val_min_ac}=={val_min_bd})")
                            print(f"      = {lhs}")
                            print(f"  RHS = I({a}=={b}) + I({c}=={d})")
                            print(f"      = {rhs}")
                            print(f"\n  Result: {lhs} < {rhs}. The condition fails for q={q}.")
                            found_counterexample = True
                            break
                    if found_counterexample: break
                if found_counterexample: break
            if found_counterexample: break

        if not found_counterexample:
            print("No counterexample found. The condition holds for a single edge.")
        print("")

# We check for q=2 and q=3.
check_fkg_condition_for_potts(3)