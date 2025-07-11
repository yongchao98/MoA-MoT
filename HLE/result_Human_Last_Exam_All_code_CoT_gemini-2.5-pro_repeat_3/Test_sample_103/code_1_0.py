def solve():
    """
    This script demonstrates a counterexample to show that f and g being
    monotone (or continuous) is not a sufficient condition for the equality
    fp(f . g) = fp(f) ∩ fp(g) to hold.
    """

    # Let L be the poset {0, 1, 2} with the usual order 0 <= 1 <= 2.
    L = {0, 1, 2}

    # Define two functions, f and g, that are monotone on L.
    # On a finite chain, monotone functions are also continuous.
    # f(x) = 1 if x<=1, 2 if x=2
    # g(x) = 0 if x<=1, 1 if x=2
    f_map = {0: 1, 1: 1, 2: 2}
    g_map = {0: 0, 1: 0, 2: 1}

    def f(x):
        return f_map[x]

    def g(x):
        return g_map[x]

    # Helper function to find the set of fixed points.
    def get_fixed_points(func, domain):
        """Returns the set of fixed points for a function over a domain."""
        return {x for x in domain if func(x) == x}

    # Helper function for composition.
    def compose(func_f, func_g):
        """Returns a new function that is the composition f(g(x))."""
        return lambda x: func_f(func_g(x))

    # Create the composed function f . g
    f_dot_g = compose(f, g)

    # Calculate the Left-Hand Side (LHS) of the equation
    fp_f_dot_g = get_fixed_points(f_dot_g, L)

    # Calculate the components of the Right-Hand Side (RHS)
    fp_f = get_fixed_points(f, L)
    fp_g = get_fixed_points(g, L)
    intersection_fp = fp_f.intersection(fp_g)

    # Print the results
    print(f"Let L = {L} be the poset.")
    print(f"f is defined as f(0)={f(0)}, f(1)={f(1)}, f(2)={f(2)}")
    print(f"g is defined as g(0)={g(0)}, g(1)={g(1)}, g(2)={g(2)}")
    print("-" * 30)
    
    # Print the final equation with the calculated sets
    print("The equation is: fp(f . g) = fp(f) ∩ fp(g)")
    print(f"LHS: fp(f . g) is the set of x where f(g(x)) = x.")
    # f(g(0)) = f(0) = 1 != 0
    # f(g(1)) = f(0) = 1 == 1  <- fixed point
    # f(g(2)) = f(1) = 1 != 2
    print(f"Calculation for LHS gives: {fp_f_dot_g}")

    print("\nRHS: fp(f) ∩ fp(g)")
    print(f"fp(f) = {fp_f}")
    print(f"fp(g) = {fp_g}")
    print(f"Calculation for RHS gives: {intersection_fp}")

    print("\n" + "-" * 30)
    # Final comparison
    print(f"Result: {fp_f_dot_g} = {intersection_fp} is {fp_f_dot_g == intersection_fp}.")
    print("This shows that f and g being monotone/continuous is not sufficient.")

solve()