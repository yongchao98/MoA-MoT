def find_counterexample_for_potts_fkg():
    """
    This script tests the FKG lattice condition for the q-state Potts model.
    A failure of this condition for a given q on a single-edge graph (deg_max=1)
    implies that the positive correlations property does not generally hold for
    any graph with maximum degree d>=1 and that number of states q.
    """
    q = 3  # Smallest number of states for which the property fails.
    s_range = range(1, q + 1)

    print(f"Checking the FKG lattice condition for the Potts model with q={q} states.")
    print("The condition must hold for all choices of states a,b,c,d from {1, ..., q}.")
    print("The condition is: I(a=b) + I(c=d) <= I(max(a,c)=max(b,d)) + I(min(a,c)=min(b,d))")
    print("-" * 20)

    for a in s_range:
        for b in s_range:
            for c in s_range:
                for d in s_range:
                    lhs = (1 if a == b else 0) + (1 if c == d else 0)
                    
                    max_ac = max(a, c)
                    max_bd = max(b, d)
                    min_ac = min(a, c)
                    min_bd = min(b, d)
                    
                    rhs = (1 if max_ac == max_bd else 0) + (1 if min_ac == min_bd else 0)
                    
                    if lhs > rhs:
                        print("Found a counterexample:")
                        print(f"Let a={a}, b={b}, c={c}, d={d}.")
                        print("\nSubstituting these values into the inequality:")
                        print(f"I({a}={b}) + I({c}={d}) <= I(max({a},{c})=max({b},{d})) + I(min({a},{c})=min({b},{d}))")
                        # Print the equation with numbers evaluated for the indicator functions
                        print(f"  {1 if a == b else 0}    +    {1 if c == d else 0}    <= I({max_ac}={max_bd}) + I({min_ac}={min_bd})")
                        print(f"      {lhs}        <= {1 if max_ac == max_bd else 0}    +    {1 if min_ac == min_bd else 0}")
                        print(f"      {lhs}        <= {rhs}")
                        print("\nThis inequality is FALSE.")
                        print("\nThis demonstrates that for q=3, the positive correlations property fails for a graph with just one edge (max degree = 1).")
                        print("Since the property must hold for ALL q>=2, it fails for any d>=1.")
                        print("The property holds for d=0 (a graph with no edges).")
                        print("Therefore, the largest integer d for which the statement is true is 0.")
                        return

    print("No counterexample found.")

# Run the demonstration
find_counterexample_for_potts_fkg()
<<<A>>>