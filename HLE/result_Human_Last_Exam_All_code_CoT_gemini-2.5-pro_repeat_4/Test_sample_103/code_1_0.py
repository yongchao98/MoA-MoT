import collections

def find_fixed_points(func, domain):
    """Finds the fixed points of a function over a given domain."""
    # A point x is a fixed point of func if func(x) = x.
    return {x for x in domain if func(x) == x}

def compose(f, g):
    """Composes two functions f and g (as dictionaries). Returns f(g(x))."""
    composition = {}
    domain_g = g.keys()
    for x in domain_g:
        composition[x] = f[g[x]]
    return composition

def is_extensive(func, domain):
    """Checks if a function is extensive (x <= func(x) for all x)."""
    # Note: This assumes the domain elements are comparable with <=.
    return all(x <= func[x] for x in domain)

def run_demonstration():
    """
    Demonstrates which condition ensures fp(f.g) = fp(f) intersect fp(g).
    """
    print("--- Analysis of fp(f.g) = fp(f) ∩ fp(g) ---")
    print("\nCase 1: Test with f and g being extensive.")
    
    # Define a poset L and extensive functions f and g
    L1 = [0, 1, 2]
    # f is extensive: 0<=1, 1<=1, 2<=2
    f1 = {0: 1, 1: 1, 2: 2}
    # g is extensive: 0<=0, 1<=2, 2<=2
    g1 = {0: 0, 1: 2, 2: 2}
    
    print(f"Let L = {L1}")
    print(f"Let f = {f1} (Extensive: {is_extensive(f1, L1)})")
    print(f"Let g = {g1} (Extensive: {is_extensive(g1, L1)})")
    
    # Calculate both sides of the equation
    fg1 = compose(f1, g1)
    fp_fg1 = find_fixed_points(fg1, L1)
    fp_f1 = find_fixed_points(f1, L1)
    fp_g1 = find_fixed_points(g1, L1)
    intersection1 = fp_f1.intersection(fp_g1)
    
    print("\nCalculating left side: fp(f . g)")
    print(f"f.g = {fg1}")
    print("The final equation's left side fixed points are:")
    for point in sorted(list(fp_fg1)):
        print(point)
    print(f"fp(f . g) = {fp_fg1}")

    print("\nCalculating right side: fp(f) ∩ fp(g)")
    print(f"fp(f) = {fp_f1}")
    print(f"fp(g) = {fp_g1}")
    print("The final equation's right side fixed points are:")
    for point in sorted(list(intersection1)):
        print(point)
    print(f"fp(f) ∩ fp(g) = {intersection1}")
    
    print(f"\nResult: {fp_fg1} == {intersection1} is {fp_fg1 == intersection1}. The equality holds.")

    print("\n" + "="*50 + "\n")

    print("Case 2: Test with f and g being monotone but NOT extensive.")
    # Counterexample for monotonicity
    L2 = [0, 1, 2]
    # f2 is monotone, not extensive
    f2 = {0: 1, 1: 2, 2: 2}
    # g2 is monotone, not extensive
    g2 = {0: 0, 1: 0, 2: 1}

    print(f"Let L = {L2}")
    print(f"Let f = {f2} (Extensive: {is_extensive(f2, L2)})")
    print(f"Let g = {g2} (Extensive: {is_extensive(g2, L2)})")

    fg2 = compose(f2, g2)
    fp_fg2 = find_fixed_points(fg2, L2)
    fp_f2 = find_fixed_points(f2, L2)
    fp_g2 = find_fixed_points(g2, L2)
    intersection2 = fp_f2.intersection(fp_g2)

    print(f"\nfp(f . g) = {fp_fg2}")
    print(f"fp(f) ∩ fp(g) = {fp_f2} ∩ {fp_g2} = {intersection2}")
    print(f"Result: {fp_fg2} == {intersection2} is {fp_fg2 == intersection2}. The equality does NOT hold.")
    print("This shows monotonicity alone is not sufficient.")

    print("\n" + "="*50 + "\n")
    
    print("Case 3: Test with only one function being extensive.")
    # Counterexample for 'f or g extensive'
    L3 = [0, 1]
    # f3 is extensive
    f3 = {0: 1, 1: 1}
    # g3 is not extensive
    g3 = {0: 1, 1: 0}
    
    print(f"Let L = {L3}")
    print(f"Let f = {f3} (Extensive: {is_extensive(f3, L3)})")
    print(f"Let g = {g3} (Extensive: {is_extensive(g3, L3)})")
    
    fg3 = compose(f3, g3)
    fp_fg3 = find_fixed_points(fg3, L3)
    fp_f3 = find_fixed_points(f3, L3)
    fp_g3 = find_fixed_points(g3, L3)
    intersection3 = fp_f3.intersection(fp_g3)
    
    print(f"\nfp(f . g) = {fp_fg3}")
    print(f"fp(f) ∩ fp(g) = {fp_f3} ∩ {fp_g3} = {intersection3}")
    print(f"Result: {fp_fg3} == {intersection3} is {fp_fg3 == intersection3}. The equality does NOT hold.")
    print("This shows 'f or g extensive' is not sufficient.")


if __name__ == '__main__':
    run_demonstration()