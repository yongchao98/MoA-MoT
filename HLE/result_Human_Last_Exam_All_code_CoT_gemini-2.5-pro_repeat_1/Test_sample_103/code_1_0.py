import collections

def demonstrate_poset_properties():
    """
    This script demonstrates the conditions under which fp(f . g) = fp(f) ∩ fp(g).
    It tests two scenarios on a simple poset to justify the answer.
    """
    # The poset L is a finite set of integers with the usual <= relation.
    L = set(range(6)) # L = {0, 1, 2, 3, 4, 5}

    def find_fp(func, domain):
        """Calculates the set of fixed points for a function on a given domain."""
        return {x for x in domain if func(x) == x}

    def print_results(case_name, f_str, g_str, f, g, f_g, domain):
        """Helper function to print the results of a test case."""
        print(f"--- {case_name} ---")
        print(f"Defining f(x) = {f_str} and g(x) = {g_str}")

        fp_f = find_fp(f, domain)
        fp_g = find_fp(g, domain)
        fp_f_g = find_fp(f_g, domain)
        intersection_fp = fp_f.intersection(fp_g)

        print(f"Poset L = {sorted(list(domain))}")
        print(f"The equation is: fp(f . g) = fp(f) ∩ fp(g)")
        print()
        print(f"fp(f) = {sorted(list(fp_f))}")
        print(f"fp(g) = {sorted(list(fp_g))}")
        print(f"RHS: fp(f) ∩ fp(g) = {sorted(list(intersection_fp))}")
        print(f"LHS: fp(f . g) = {sorted(list(fp_f_g))}")

        if fp_f_g == intersection_fp:
            print("\nResult: The equality holds.")
        else:
            print("\nResult: The equality DOES NOT hold.")
        print("-" * (len(case_name) + 8))
        print()


    # Case 1: Test the condition "f and g are monotone"
    # This serves as a counterexample for options A, C, F, G.
    def g_mono(n):
        return min(n + 1, max(L))
    def f_mono(n):
        return max(n - 1, min(L))
    def f_dot_g_mono(n):
        return f_mono(g_mono(n))

    print_results("Case 1: Monotone Functions (Counterexample)",
                  "max(x - 1, 0)", "min(x + 1, 5)",
                  f_mono, g_mono, f_dot_g_mono, L)

    # Case 2: Test the condition "f and g are extensive"
    # An extensive function h satisfies x <= h(x) for all x.
    def g_extensive(n):
        return min(n + 2, max(L))
    def f_extensive(n):
        return min(n + 1, max(L))
    def f_dot_g_extensive(n):
        return f_extensive(g_extensive(n))

    print_results("Case 2: Extensive Functions",
                  "min(x + 1, 5)", "min(x + 2, 5)",
                  f_extensive, g_extensive, f_dot_g_extensive, L)

if __name__ == '__main__':
    demonstrate_poset_properties()