def get_fp(func, domain):
    """Computes the set of fixed points for a function."""
    return {x for x in domain if func(x) == x}

def compose(f, g):
    """Composes two functions represented as dictionaries."""
    return lambda x: f[g[x]]

def run_demonstration():
    """
    Demonstrates the conditions for the fixed point equality fp(f.g) = fp(f) ∩ fp(g).
    """
    # The poset L is the set {0, 1, 2} with the usual order 0 <= 1 <= 2.
    L = {0, 1, 2}
    print("This script demonstrates the conditions for fp(f . g) = fp(f) ∩ fp(g).")
    print(f"We use the poset L = {L} with the usual order 0 <= 1 <= 2.\n")

    print("--- Case 1: Counterexample (g is extensive, f is not) ---")
    
    # Define functions as dictionaries for the counterexample
    # g1 is extensive (x <= g1(x) for all x)
    # f1 is not extensive because f1(1) < 1
    f1 = {0: 0, 1: 0, 2: 2}
    g1 = {0: 1, 1: 1, 2: 2}
    
    print(f"Let f = {f1}")
    print(f"Let g = {g1}\n")

    # Calculate fixed points and their intersection (RHS of the equation)
    fp_f1 = get_fp(lambda x: f1[x], L)
    fp_g1 = get_fp(lambda x: g1[x], L)
    intersection_fp = fp_f1.intersection(fp_g1)
    
    # Calculate fixed points of the composition (LHS of the equation)
    f1_g1_comp = compose(f1, g1)
    fp_comp = get_fp(f1_g1_comp, L)

    print("The equation to check is: fp(f . g) = fp(f) ∩ fp(g)")
    print(f"LHS: fp(f . g) = {sorted(list(fp_comp))}")
    print(f"RHS: fp(f) ∩ fp(g) = {sorted(list(intersection_fp))}")
    print(f"Substituting the values: {sorted(list(fp_comp))} = {sorted(list(intersection_fp))}, which is {fp_comp == intersection_fp}.")
    print("This demonstrates that 'f or g extensive' is not a sufficient condition.\n")
    print("-" * 50)

    print("\n--- Case 2: Sufficient Condition (f and g are extensive) ---")
    
    # Define two extensive functions
    # f2 is extensive: 0<=0, 1<=2, 2<=2
    # g2 is extensive: 0<=1, 1<=1, 2<=2
    f2 = {0: 0, 1: 2, 2: 2}
    g2 = {0: 1, 1: 1, 2: 2}
    
    print(f"Let f = {f2}")
    print(f"Let g = {g2}\n")

    # Calculate fixed points and their intersection (RHS)
    fp_f2 = get_fp(lambda x: f2[x], L)
    fp_g2 = get_fp(lambda x: g2[x], L)
    intersection_fp_2 = fp_f2.intersection(fp_g2)
    
    # Calculate fixed points of the composition (LHS)
    f2_g2_comp = compose(f2, g2)
    fp_comp_2 = get_fp(f2_g2_comp, L)

    print("The equation to check is: fp(f . g) = fp(f) ∩ fp(g)")
    print(f"LHS: fp(f . g) = {sorted(list(fp_comp_2))}")
    print(f"RHS: fp(f) ∩ fp(g) = {sorted(list(intersection_fp_2))}")
    print(f"Substituting the values: {sorted(list(fp_comp_2))} = {sorted(list(intersection_fp_2))}, which is {fp_comp_2 == intersection_fp_2}.")
    print("This demonstrates that when f and g are both extensive, the equality holds.")

if __name__ == "__main__":
    run_demonstration()