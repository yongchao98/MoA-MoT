def check_property():
    """
    This script demonstrates a counterexample for the statement:
    If f and g are monotone, then fp(f . g) = fp(f) ∩ fp(g).
    """
    # The poset is the set L = {0, 1, 2} with the usual order 0 <= 1 <= 2.
    L = [0, 1, 2]

    # Define two monotone functions f and g on L.
    # We use dictionaries to represent the functions.
    # g(0)=1, g(1)=2, g(2)=2. This is monotone.
    g = {0: 1, 1: 2, 2: 2}
    # f(0)=0, f(1)=0, f(2)=1. This is monotone.
    f = {0: 0, 1: 0, 2: 1}

    print(f"Poset L = {L}")
    print(f"Function g = {g}")
    print(f"Function f = {f}")
    print("-" * 20)

    # 1. Calculate fp(f) ∩ fp(g)
    fp_f = {x for x in L if f[x] == x}
    fp_g = {x for x in L if g[x] == x}
    intersection_fp = fp_f.intersection(fp_g)

    print(f"fp(f) = {{{', '.join(map(str, sorted(list(fp_f))))}}}")
    print(f"fp(g) = {{{', '.join(map(str, sorted(list(fp_g))))}}}")
    print(f"fp(f) ∩ fp(g) = {{{', '.join(map(str, sorted(list(intersection_fp))))}}}")
    print("-" * 20)

    # 2. Calculate fp(f . g)
    # (f . g)(x) = f(g(x))
    fp_f_g = {x for x in L if f[g[x]] == x}

    print(f"Calculating fixed points of f(g(x)):")
    for x in L:
        print(f"f(g({x})) = f({g[x]}) = {f[g[x]]}")

    print(f"fp(f . g) = {{{', '.join(map(str, sorted(list(fp_f_g))))}}}")
    print("-" * 20)

    # 3. Compare the two sets
    are_equal = (fp_f_g == intersection_fp)
    print(f"Is fp(f . g) == fp(f) ∩ fp(g)? {are_equal}")
    if not are_equal:
        print("The equality does not hold. 'f and g monotone' is not a sufficient condition.")

check_property()