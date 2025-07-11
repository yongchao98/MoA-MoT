def analyze_fixed_point_conditions():
    """
    Analyzes various conditions on functions f and g to determine when
    the equality fp(f . g) = fp(f) ∩ fp(g) holds for a poset (L, ≤).
    
    The code tests several hypotheses by attempting to build counterexamples.
    fp(f) denotes the set of fixed points of f.
    f . g is function composition f(g(x)).
    """
    
    # --- Helper functions ---
    
    def fp(f, domain):
        """Computes the set of fixed points of a function."""
        return {x for x in domain if f(x) == x}

    def compose(f, g):
        """Returns the composite function f(g(x))."""
        return lambda x: f(g(x))

    def is_extensive(f, domain, leq):
        """Checks if a function is extensive (x ≤ f(x) for all x)."""
        for x in domain:
            if not leq(x, f(x)):
                return False
        return True

    print("Investigating the condition for fp(f.g) = fp(f) ∩ fp(g)\n")
    print("We will test several options by constructing counterexamples.")
    print("-" * 50)

    # --- Test Case 1: Monotonicity (Counterexample) ---
    print("Hypothesis: 'f and g are monotone' is sufficient.")
    L1 = {1, 2, 3}
    leq1 = lambda x, y: x <= y  # Standard integer order
    f1_dict = {1: 1, 2: 1, 3: 2}
    g1_dict = {1: 2, 2: 3, 3: 3}
    f1 = lambda x: f1_dict[x]
    g1 = lambda x: g1_dict[x]

    fp_f1 = fp(f1, L1)
    fp_g1 = fp(g1, L1)
    intersection1 = fp_f1.intersection(fp_g1)
    fg1 = compose(f1, g1)
    fp_fg1 = fp(fg1, L1)

    print(f"On poset L={L1} with f={f1_dict}, g={g1_dict}:")
    print(f"Equation: fp(f.g) == fp(f) ∩ fp(g)")
    print(f"fp(f.g) is: {fp_fg1}")
    print(f"fp(f) ∩ fp(g) is: {fp_f1} ∩ {fp_g1} = {intersection1}")
    print(f"Are they equal? {'Yes' if fp_fg1 == intersection1 else 'No'}")
    print("Result: Monotonicity is not sufficient.\n")
    print("-" * 50)

    # --- Test Case 2: Only one function is extensive (Counterexample) ---
    print("Hypothesis: 'f or g is extensive' is sufficient.")
    L2 = {1, 2}
    leq2 = lambda x, y: x <= y
    
    # Case 2a: Only f is extensive
    f2a_dict = {1: 2, 2: 2} # Extensive
    g2a_dict = {1: 1, 2: 1} # Not extensive
    f2a = lambda x: f2a_dict[x]
    g2a = lambda x: g2a_dict[x]
    fp_f2a = fp(f2a, L2)
    fp_g2a = fp(g2a, L2)
    intersection2a = fp_f2a.intersection(fp_g2a)
    fg2a = compose(f2a, g2a)
    fp_fg2a = fp(fg2a, L2)

    print(f"Case 'f is extensive, g is not' on L={L2} with f={f2a_dict}, g={g2a_dict}:")
    print(f"Equation: fp(f.g) == fp(f) ∩ fp(g)")
    print(f"fp(f.g) is: {fp_fg2a}")
    print(f"fp(f) ∩ fp(g) is: {fp_f2a} ∩ {fp_g2a} = {intersection2a}")
    print(f"Are they equal? {'Yes' if fp_fg2a == intersection2a else 'No'}")
    print("Result: 'f or g extensive' is not sufficient.\n")
    print("-" * 50)
    
    # --- Test Case 3: Both functions are extensive (Illustrative Example + Proof) ---
    print("Hypothesis: 'f and g are extensive' is sufficient.")
    L3 = {0, 1, 2, 3}
    leq3 = lambda x, y: x <= y
    f3_dict = {0: 1, 1: 1, 2: 3, 3: 3} # Extensive: x <= f(x)
    g3_dict = {0: 0, 1: 2, 2: 2, 3: 3} # Extensive: x <= g(x)
    f3 = lambda x: f3_dict[x]
    g3 = lambda x: g3_dict[x]
    
    print(f"On poset L={L3} with extensive f={f3_dict}, g={g3_dict}:")
    fp_f3 = fp(f3, L3)
    fp_g3 = fp(g3, L3)
    intersection3 = fp_f3.intersection(fp_g3)
    fg3 = compose(f3, g3)
    fp_fg3 = fp(fg3, L3)
    
    print(f"Equation: fp(f.g) == fp(f) ∩ fp(g)")
    print(f"fp(f.g) is: {fp_fg3}")
    print(f"fp(f) ∩ fp(g) is: {fp_f3} ∩ {fp_g3} = {intersection3}")
    print(f"Are they equal? {'Yes' if fp_fg3 == intersection3 else 'No'}")

    print("\nFormal Proof Sketch:")
    print("Assume f and g are extensive (x ≤ f(x) and x ≤ g(x) for all x).")
    print("Let x be in fp(f.g), which means f(g(x)) = x.")
    print("1. By g's extensivity: x ≤ g(x).")
    print("2. Let y = g(x). Then x ≤ y and f(y) = x.")
    print("3. By f's extensivity: y ≤ f(y).")
    print("4. Substituting from step 2: g(x) ≤ x.")
    print("5. From steps 1 and 4, and antisymmetry of ≤, we have g(x) = x.")
    print("6. Substitute g(x)=x into f(g(x))=x to get f(x) = x.")
    print("So, x is a fixed point of both f and g. This proves fp(f.g) ⊆ fp(f) ∩ fp(g).")
    print("Result: This condition is sufficient.")

analyze_fixed_point_conditions()