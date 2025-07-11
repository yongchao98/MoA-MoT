def get_fixed_points(func, domain):
    """Finds the fixed points of a function over a given domain."""
    return {x for x in domain if func(x) == x}

def main():
    """
    Demonstrates that for two extensive functions f and g,
    fp(f . g) = fp(f) ∩ fp(g).
    """
    # Let L be the poset {0, 1, 2, 3, 4, 5} with the usual order <=
    L = set(range(6))

    # Define an extensive function f. f(x) >= x for all x in L.
    # f(x) = x + 1 for x < 5, and f(5) = 5.
    def f(x):
        return min(x + 1, 5)

    # Define another extensive function g. g(x) >= x for all x in L.
    # g(0)=0, g(1)=2, g(2)=4, g(3)=4, g(4)=4, g(5)=5
    def g(x):
        if x == 0: return 0
        if x == 1: return 2
        if x in [2,3,4]: return 4
        if x == 5: return 5
        
    # Check for extensivity (optional, for verification)
    is_f_extensive = all(x <= f(x) for x in L)
    is_g_extensive = all(x <= g(x) for x in L)
    print(f"Is f extensive? {is_f_extensive}")
    print(f"Is g extensive? {is_g_extensive}\n")
    
    # Define the composition f . g
    def f_compose_g(x):
        return f(g(x))

    # Calculate the set of fixed points for each function
    fp_f = get_fixed_points(f, L)
    fp_g = get_fixed_points(g, L)
    fp_f_compose_g = get_fixed_points(f_compose_g, L)

    # Calculate the intersection of fp(f) and fp(g)
    intersection_fp = fp_f.intersection(fp_g)
    
    # The final equation is fp(f.g) = fp(f) ∩ fp(g)
    # We will print each component of the equation
    print("Evaluating the equation: fp(f.g) = fp(f) ∩ fp(g)")
    print("-" * 50)
    print(f"fp(f.g) is the set: {fp_f_compose_g}")
    print(f"fp(f) is the set: {fp_f}")
    print(f"fp(g) is the set: {fp_g}")
    print(f"fp(f) ∩ fp(g) is the set: {intersection_fp}")
    print("-" * 50)

    # Verify and print the result
    are_equal = (fp_f_compose_g == intersection_fp)
    print(f"Is fp(f.g) equal to fp(f) ∩ fp(g)? {are_equal}")

if __name__ == "__main__":
    main()