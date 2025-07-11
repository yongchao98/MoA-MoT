def check_fixed_point_equality():
    """
    Demonstrates the condition for fp(f.g) = fp(f) intersect fp(g).
    """
    # Let L be a poset, represented by a range of integers with the usual <= order.
    L = range(20)

    # --- Case 1: f and g are extensive ---
    # An extensive function h satisfies x <= h(x).
    # f(x) = x + 2 (but capped at 19 to stay within L)
    # g(x) = x + 1 (but capped at 19 to stay within L)
    f = lambda x: min(x + 2, 19)
    g = lambda x: min(x + 1, 19)
    
    # Check if f and g are extensive on L
    is_f_extensive = all(x <= f(x) for x in L)
    is_g_extensive = all(x <= g(x) for x in L)

    print("--- Case 1: f and g are extensive ---")
    print(f"Is f(x) = min(x+2, 19) extensive? {is_f_extensive}")
    print(f"Is g(x) = min(x+1, 19) extensive? {is_g_extensive}")

    # Compute fixed points
    fp_f = {x for x in L if f(x) == x}
    fp_g = {x for x in L if g(x) == x}
    fp_fg = {x for x in L if f(g(x)) == x}
    
    # Compute the intersection
    intersection = fp_f.intersection(fp_g)

    print("\nEquation: fp(f . g) = fp(f) \u2229 fp(g)")
    print(f"fp(f . g) = {fp_fg}")
    print(f"fp(f) \u2229 fp(g) = {fp_f} \u2229 {fp_g} = {intersection}")
    print(f"Does the equality hold? {fp_fg == intersection}")
    
    # --- Case 2: One function is not extensive ---
    # h(x) = x - 5 (but capped at 0) is not extensive for x > 0.
    h = lambda x: max(x - 5, 0)
    
    is_h_extensive = all(x <= h(x) for x in L)
    
    print("\n\n--- Case 2: g is replaced by non-extensive function h ---")
    print(f"Is h(x) = max(x-5, 0) extensive? {is_h_extensive}")
    
    # Compute fixed points for f and h
    fp_h = {x for x in L if h(x) == x}
    fp_fh = {x for x in L if f(h(x)) == x}

    # Compute the intersection
    intersection_fh = fp_f.intersection(fp_h)
    
    print("\nEquation: fp(f . h) = fp(f) \u2229 fp(h)")
    print(f"fp(f . h) = {fp_fh}")
    print(f"fp(f) \u2229 fp(h) = {fp_f} \u2229 {fp_h} = {intersection_fh}")
    print(f"Does the equality hold? {fp_fh == intersection_fh}")

check_fixed_point_equality()