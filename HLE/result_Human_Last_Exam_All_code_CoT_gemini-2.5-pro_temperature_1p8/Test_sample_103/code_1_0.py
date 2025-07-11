def check_property():
    """
    This script demonstrates that "f and g are monotone" is not a sufficient
    condition for fp(f*g) = fp(f) ∩ fp(g).
    It sets up a counterexample on the poset L = {0, 1} with 0 <= 1.
    """
    
    # Define the poset L
    L = {0, 1}

    # Define two functions f and g that are monotone but not extensive.
    # f(0)=1, f(1)=1
    # g(0)=0, g(1)=0
    # Both are monotone because for x <= y, f(x) <= f(y) and g(x) <= g(y).
    # e.g., 0 <= 1, f(0)=1 <= f(1)=1.
    f_map = {0: 1, 1: 1}
    g_map = {0: 0, 1: 0}
    
    # Define a helper to calculate fixed points
    def get_fixed_points(func_map, domain):
        return {x for x in domain if func_map[x] == x}

    # Calculate fixed points for f and g
    fp_f = get_fixed_points(f_map, L)
    fp_g = get_fixed_points(g_map, L)

    # Calculate the intersection of the fixed point sets
    fp_intersection = fp_f.intersection(fp_g)

    # Define the composition f*g (f(g(x)))
    f_circ_g_map = {x: f_map[g_map[x]] for x in L}

    # Calculate fixed points of the composed function
    fp_f_circ_g = get_fixed_points(f_circ_g_map, L)

    print("This script checks if fp(f*g) == fp(f) ∩ fp(g) for two monotone functions.")
    print(f"Let L = {L} with 0 <= 1.")
    print(f"Let f be defined as f(x) = {f_map}.")
    print(f"Let g be defined as g(x) = {g_map}.")
    print("\n--- Equation ---")
    print("fp(f * g) = fp(f) ∩ fp(g)\n")
    
    print("--- Left Hand Side ---")
    # Using sorted list for consistent output order
    print(f"The fixed points of the composition f*g are: {sorted(list(fp_f_circ_g))}")
    print("\n--- Right Hand Side ---")
    print(f"The fixed points of f are: {sorted(list(fp_f))}")
    print(f"The fixed points of g are: {sorted(list(fp_g))}")
    # Using sorted list for consistent output order
    intersection_list = sorted(list(fp_intersection)) if fp_intersection else "the empty set"
    print(f"The intersection of their fixed points is: {intersection_list}\n")

    if fp_f_circ_g == fp_intersection:
        print("Conclusion: The equality holds for this example.")
    else:
        print("Conclusion: The equality does NOT hold. Monotonicity is not a sufficient condition.")

check_property()