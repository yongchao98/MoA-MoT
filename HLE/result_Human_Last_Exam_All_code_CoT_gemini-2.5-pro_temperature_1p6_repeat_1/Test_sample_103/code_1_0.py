def find_fixed_points(func, domain):
    """Finds the fixed points of a function over a given domain."""
    return {x for x in domain if func(x) == x}

def test_case(description, domain, f, g):
    """
    Tests the equality fp(f*g) = fp(f) ∩ fp(g) for a given scenario.
    """
    print(f"--- {description} ---")

    # The composite function f(g(x))
    f_g = lambda x: f(g(x))

    # Calculate the fixed point sets
    fp_f = find_fixed_points(f, domain)
    fp_g = find_fixed_points(g, domain)
    fp_f_g = find_fixed_points(f_g, domain)

    # Calculate the intersection
    intersection = fp_f.intersection(fp_g)

    # Print the sets
    print(f"Domain L = {set(domain)}")
    print(f"fp(f) = {fp_f}")
    print(f"fp(g) = {fp_g}")
    print(f"fp(f) ∩ fp(g) = {intersection}")
    print(f"fp(f * g) = {fp_f_g}\n")
    
    # Print the equation with the computed values
    print("Checking the equation: fp(f * g) = fp(f) ∩ fp(g)")
    print(f"Result: {fp_f_g} = {intersection}")

    # Check for equality
    if fp_f_g == intersection:
        print("The equality holds true.")
    else:
        print("The equality does not hold.")
    print("-" * (len(description) + 6) + "\n")

def main():
    # Define the poset (a set of integers with the usual <= order)
    L = list(range(5))

    # Case 1: Both f and g are extensive.
    # An extensive function h(x) satisfies x <= h(x).
    f1 = lambda x: min(x + 1, 4) # f1 is extensive
    g1 = lambda x: min(x + 2, 4) # g1 is extensive
    test_case("Case 1: f and g are both extensive", L, f1, g1)

    # Case 2: A counterexample (e.g., for 'f or g extensive')
    # Here g2 is extensive, but f2 is not. This also serves as a
    # counterexample for 'f and g monotone' since both are monotone.
    f2 = lambda x: max(x - 1, 0) # f2 is monotone, but not extensive
    g2 = lambda x: min(x + 1, 4) # g2 is monotone and extensive
    test_case("Case 2: Counterexample (g is extensive, f is not)", L, f2, g2)


if __name__ == "__main__":
    main()