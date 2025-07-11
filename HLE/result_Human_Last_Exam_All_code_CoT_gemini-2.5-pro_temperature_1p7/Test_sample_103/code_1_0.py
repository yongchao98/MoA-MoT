def compute_fp(func, poset):
    """Computes the set of fixed points for a function on a poset."""
    return {x for x in poset if func(x) == x}

def compose(f, g):
    """Returns the composition of two functions."""
    return lambda x: f(g(x))

def run_demonstration():
    """
    Demonstrates that 'f and g monotone' is not a sufficient condition.
    This serves as a counterexample for several weaker options.
    """
    # Let L be a poset {0, 1, 2, 3} with the standard order <=
    L = {0, 1, 2, 3}
    print(f"Poset L = {L}")
    print("-" * 20)
    
    # Define two monotone functions f and g on L.
    # Note: g is also extensive, while f is not.
    f_map = {0: 0, 1: 1, 2: 1, 3: 3}
    g_map = {0: 0, 1: 2, 2: 3, 3: 3}
    
    f = lambda x: f_map[x]
    g = lambda x: g_map[x]

    # f is monotone because for x <= y, f(x) <= f(y) (0<=1, 1<=1, 1<=3)
    # g is monotone because for x <= y, g(x) <= g(y) (0<=2, 2<=3, 3<=3)
    print("f and g are constructed to be monotone.")
    print(f"f mapping: {f_map}")
    print(f"g mapping: {g_map}")
    print("-" * 20)

    # Calculate fp(f) ∩ fp(g)
    fp_f = compute_fp(f, L)
    fp_g = compute_fp(g, L)
    intersection_of_fps = fp_f.intersection(fp_g)

    print(f"fp(f) = {fp_f}")
    print(f"fp(g) = {fp_g}")
    print(f"The equation side 'fp(f) ∩ fp(g)' evaluates to: {intersection_of_fps}")
    print("-" * 20)
    
    # Calculate fp(f . g)
    f_dot_g = compose(f, g)
    fp_composition = compute_fp(f_dot_g, L)
    
    # (f.g)(0) = f(g(0)) = f(0) = 0.  (Fixed point)
    # (f.g)(1) = f(g(1)) = f(2) = 1.  (Fixed point)
    # (f.g)(2) = f(g(2)) = f(3) = 3.  (Not a fixed point)
    # (f.g)(3) = f(g(3)) = f(3) = 3.  (Fixed point)
    
    print(f"The function '(f . g)' has fixed points: {fp_composition}")
    print(f"The equation side 'fp(f . g)' evaluates to: {fp_composition}")
    print("-" * 20)

    # Compare the two sets
    are_equal = (fp_composition == intersection_of_fps)
    print(f"Result of fp(f . g) = fp(f) ∩ fp(g): {are_equal}")
    if not are_equal:
        print("Since the sets are not equal, the condition 'f and g monotone' is not sufficient.")

if __name__ == '__main__':
    run_demonstration()
