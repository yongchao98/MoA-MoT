def main():
    """
    This script demonstrates that "f and g monotone" is not a sufficient
    condition for fp(f . g) = fp(f) ∩ fp(g).
    """
    # Let L be a finite subset of integers, a poset with the usual <= order.
    L = list(range(-5, 6))

    # Let f and g be two monotone functions.
    # f is not extensive, g is extensive.
    def f(x):
        return x - 1

    def g(x):
        return x + 1

    # Calculate fp(f . g) = {x in L | f(g(x)) = x}
    fp_f_g = []
    for x in L:
        # We need to ensure g(x) is in L for f(g(x)) to be well-defined within our set.
        if g(x) in L and f(g(x)) == x:
            fp_f_g.append(x)

    # Calculate fp(f) = {x in L | f(x) = x}
    fp_f = []
    for x in L:
        if f(x) == x:
            fp_f.append(x)

    # Calculate fp(g) = {x in L | g(x) = x}
    fp_g = []
    for x in L:
        if g(x) == x:
            fp_g.append(x)

    # Calculate fp(f) ∩ fp(g)
    set_fp_g = set(fp_g)
    intersection = [val for val in fp_f if val in set_fp_g]

    # Print the results
    # For f(g(x)) = f(x+1) = x, all x where x and x+1 are in L will be fixed points.
    # In L = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5], g(x) is in L for x from -5 to 4.
    expected_fp_f_g = list(range(-5, 5))
    
    print(f"Let f(x) = x - 1 and g(x) = x + 1.")
    print(f"Over the domain L = {L}:")
    print(f"fp(f . g) is the set of x where f(g(x)) = x.")
    print(f"f(g(x)) = (x + 1) - 1 = x.")
    print(f"fp(f . g) = {expected_fp_f_g}")
    print(f"fp(f) is the set of x where f(x) = x (i.e., x - 1 = x).")
    print(f"fp(f) = {fp_f}")
    print(f"fp(g) is the set of x where g(x) = x (i.e., x + 1 = x).")
    print(f"fp(g) = {fp_g}")
    print(f"fp(f) ∩ fp(g) = {intersection}")
    
    print("\nFinal Equation Check:")
    print(f"{expected_fp_f_g} = {intersection}  <-- This is FALSE.")
    print("Therefore, the condition is not sufficient.")


if __name__ == "__main__":
    main()
