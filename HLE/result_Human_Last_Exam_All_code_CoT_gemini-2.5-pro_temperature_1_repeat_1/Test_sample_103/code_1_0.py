import functools

def get_fp(func, domain):
    """Computes the set of fixed points for a function."""
    return {x for x in domain if func(x) == x}

def is_extensive(func, domain, leq):
    """Checks if a function is extensive (x <= f(x) for all x)."""
    return all(leq(x, func(x)) for x in domain)

def main():
    """
    Demonstrates the condition for fp(f.g) = fp(f) intersect fp(g).
    """
    # Define the poset (L, <=)
    # L is the set of elements, and the order is the standard integer comparison.
    L = list(range(4))  # L = {0, 1, 2, 3}
    leq = lambda a, b: a <= b

    print("--- Scenario 1: f and g are both extensive ---")
    # Define two extensive functions
    # f(x) = x | 2 (bitwise OR) is extensive on {0,1,2,3}
    # g(x) = x | 1 (bitwise OR) is extensive on {0,1,2,3}
    f = lambda x: x | 2
    g = lambda x: x | 1
    
    # Verify properties
    print(f"Is f extensive? {is_extensive(f, L, leq)}")
    print(f"Is g extensive? {is_extensive(g, L, leq)}\n")

    # Compute the sets for the equality
    fp_f = get_fp(f, L)
    fp_g = get_fp(g, L)
    intersection = fp_f.intersection(fp_g)
    
    f_dot_g = lambda x: f(g(x))
    fp_f_dot_g = get_fp(f_dot_g, L)

    # Print the equation and the result
    print("Equation: fp(f . g) = fp(f) ∩ fp(g)")
    print(f"fp(f . g) = {sorted(list(fp_f_dot_g))}")
    print(f"fp(f) = {sorted(list(fp_f))}, fp(g) = {sorted(list(fp_g))}")
    print(f"fp(f) ∩ fp(g) = {sorted(list(intersection))}")
    print(f"Result: The equality holds. ({sorted(list(fp_f_dot_g))} == {sorted(list(intersection))})")

    print("\n" + "="*50 + "\n")

    print("--- Scenario 2: f is extensive, g is not ---")
    # Define f as extensive, g as not extensive
    # f(x) = x | 2 is extensive
    # g(x) = x & 1 (bitwise AND) is NOT extensive (e.g., 2 <= g(2)=0 is False)
    f2 = lambda x: x | 2
    g2 = lambda x: x & 1

    # Verify properties
    print(f"Is f2 extensive? {is_extensive(f2, L, leq)}")
    print(f"Is g2 extensive? {is_extensive(g2, L, leq)}\n")
    
    # Compute the sets for the equality
    fp_f2 = get_fp(f2, L)
    fp_g2 = get_fp(g2, L)
    intersection2 = fp_f2.intersection(fp_g2)
    
    f2_dot_g2 = lambda x: f2(g2(x))
    fp_f2_dot_g2 = get_fp(f2_dot_g2, L)

    # Print the equation and the result
    print("Equation: fp(f . g) = fp(f) ∩ fp(g)")
    print(f"fp(f . g) = {sorted(list(fp_f2_dot_g2))}")
    print(f"fp(f) = {sorted(list(fp_f2))}, fp(g) = {sorted(list(fp_g2))}")
    print(f"fp(f) ∩ fp(g) = {sorted(list(intersection2))}")
    print(f"Result: The equality fails. ({sorted(list(fp_f2_dot_g2))} != {sorted(list(intersection2))})")

if __name__ == "__main__":
    main()