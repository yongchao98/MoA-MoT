import itertools

def check_holley_inequality(q):
    """
    Checks if the Holley criterion inequality holds for the q-state Potts model.

    The inequality is:
    I(max(a,b)==max(c,d)) + I(min(a,b)==min(c,d)) >= I(a==c) + I(b==d)
    where a, b, c, d are spin values in {1, ..., q}.
    """
    print(f"--- Checking for q = {q} ---")
    spins = range(1, q + 1)
    # Iterate through all combinations of 4 spins
    for a, b, c, d in itertools.product(spins, repeat=4):
        # Left Hand Side (LHS) of the inequality
        lhs = (1 if max(a, b) == max(c, d) else 0) + \
              (1 if min(a, b) == min(c, d) else 0)

        # Right Hand Side (RHS) of the inequality
        rhs = (1 if a == c else 0) + \
              (1 if b == d else 0)

        if lhs < rhs:
            print("Holley condition FAILS.")
            print(f"Counterexample found for (a,b,c,d) = ({a}, {b}, {c}, {d})")
            print(f"LHS = I(max({a},{b})==max({c},{d})) + I(min({a},{b})==min({c},{d}))")
            print(f"    = I({max(a,b)}=={max(c,d)}) + I({min(a,b)}=={min(c,d)})")
            print(f"    = {1 if max(a,b)==max(c,d) else 0} + {1 if min(a,b)==min(c,d) else 0} = {lhs}")
            print(f"RHS = I({a}=={c}) + I({b}=={d})")
            print(f"    = {1 if a==c else 0} + {1 if b==d else 0} = {rhs}")
            print(f"Inequality Fails: {lhs} >= {rhs} is FALSE.")
            return False

    print("Holley condition HOLDS for all spin combinations.")
    return True

# Check for q=2
check_holley_inequality(2)

print("\n")

# Check for q=3
check_holley_inequality(3)

print("\nAs shown, the inequality required for the positive correlations property fails for q=3.")
print("This failure occurs on any single edge, meaning any graph with a maximum degree d>=1 will not satisfy the property for all q.")
print("For the property to hold for all q>=2, the graph must have no edges. The only connected graph with no edges is the single-vertex graph, which has maximum degree 0.")
print("Therefore, the largest value for d is 0.")
