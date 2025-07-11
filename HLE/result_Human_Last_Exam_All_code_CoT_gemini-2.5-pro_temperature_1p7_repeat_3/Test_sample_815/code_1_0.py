import math

def get_order_GL(n, q):
    """Calculates the order of the General Linear Group GL(n, q)."""
    if q == 1:
        return 1
    order = 1
    for i in range(n):
        order *= (q**n - q**i)
    return order

def get_order_SL(n, q):
    """Calculates the order of the Special Linear Group SL(n, q)."""
    if q == 1:
        return 1
    return get_order_GL(n, q) // (q - 1)

def calculate_involutions():
    """
    Calculates and compares the number of involutions for PSL(3,9) and PSL(4,3).
    """
    # --- Calculation for PSL(3,9) ---
    n1, q1 = 3, 9
    print("Calculating the number of involutions for PSL(3,9):")
    # The center of SL(3,9) is trivial since gcd(3, 9-1) = gcd(3,8) = 1.
    # Therefore, PSL(3,9) is the same as SL(3,9).
    # We need to count the elements of order 2 in SL(3,9).
    # In odd characteristic, for n=3, involutions must have eigenvalues (1, -1, -1) to have determinant 1.
    # These form a single conjugacy class. Its size is |SL(3,9)| / |C|, where C is the centralizer.
    
    order_sl_3_9 = get_order_SL(n1, q1)
    
    # The centralizer in GL(3,9) is GL(1,9) x GL(2,9).
    order_gl_1_9 = get_order_GL(1, q1)
    order_gl_2_9 = get_order_GL(2, q1)
    centralizer_gl_size_1 = order_gl_1_9 * order_gl_2_9
    
    # The centralizer in SL(3,9) has size |C_GL| / (q-1).
    centralizer_sl_size_1 = centralizer_gl_size_1 // (q1 - 1)
    
    inv_psl_3_9 = order_sl_3_9 // centralizer_sl_size_1
    
    print("The number of involutions is given by the size of a conjugacy class in SL(3,9):")
    print(f"|SL(3,9)| = {order_sl_3_9}")
    print(f"Centralizer size in GL(3,9) = |GL(1,9)| * |GL(2,9)| = {order_gl_1_9} * {order_gl_2_9} = {centralizer_gl_size_1}")
    print(f"Centralizer size in SL(3,9) = {centralizer_gl_size_1} / (9 - 1) = {centralizer_sl_size_1}")
    print(f"Number of involutions in PSL(3,9) = {order_sl_3_9} / {centralizer_sl_size_1} = {inv_psl_3_9}")
    print("-" * 30)
    
    # --- Calculation for PSL(4,3) ---
    n2, q2 = 4, 3
    print("Calculating the number of involutions for PSL(4,3):")
    
    # The center Z(SL(4,3)) has size gcd(4, 3-1) = 2. Let Z = {I, -I}.
    # Preimages g of involutions in PSL(4,3) satisfy g^2 is in Z and g is not in Z.
    # Case 1: g^2 = I (g is a non-central involution in SL(4,3)).
    #         These have eigenvalues (1,1,-1,-1) to have det 1.
    order_sl_4_3 = get_order_SL(n2, q2)
    order_gl_2_3 = get_order_GL(2, q2)
    centralizer_gl_size_2a = order_gl_2_3 * order_gl_2_3
    centralizer_sl_size_2a = centralizer_gl_size_2a // (q2 - 1)
    num_g_sq_I = order_sl_4_3 // centralizer_sl_size_2a
    
    print("First, we count preimages g in SL(4,3).")
    print("Case 1: g is not central and g^2 = I.")
    print("Number of such elements = |SL(4,3)| / |C_SL(t)| for t=diag(1,1,-1,-1).")
    print(f"|SL(4,3)| = {order_sl_4_3}")
    print(f"|C_GL(t)| = |GL(2,3)| * |GL(2,3)| = {order_gl_2_3} * {order_gl_2_3} = {centralizer_gl_size_2a}")
    print(f"|C_SL(t)| = {centralizer_gl_size_2a} / (3 - 1) = {centralizer_sl_size_2a}")
    print(f"Number of preimages with g^2 = I is {order_sl_4_3} / {centralizer_sl_size_2a} = {num_g_sq_I}")
    
    # Case 2: g^2 = -I.
    #         Centralizer in GL(4,3) is isomorphic to GL(2, 3^2) = GL(2,9).
    centralizer_gl_size_2b = order_gl_2_9
    centralizer_sl_size_2b = centralizer_gl_size_2b // (q2 - 1)
    num_g_sq_neg_I = order_sl_4_3 // centralizer_sl_size_2b
    
    print("\nCase 2: g^2 = -I.")
    print("Number of such elements = |SL(4,3)| / |C_SL(g)|.")
    print(f"Centralizer C_GL(g) is GL(2,9), with order {centralizer_gl_size_2b}")
    print(f"|C_SL(g)| = {centralizer_gl_size_2b} / (3 - 1) = {centralizer_sl_size_2b}")
    print(f"Number of preimages with g^2 = -I is {order_sl_4_3} / {centralizer_sl_size_2b} = {num_g_sq_neg_I}")
    
    # Total involutions in PSL(4,3)
    total_preimages = num_g_sq_I + num_g_sq_neg_I
    center_size = math.gcd(n2, q2 - 1)
    inv_psl_4_3 = total_preimages // center_size
    
    print(f"\nTotal preimages = {num_g_sq_I} + {num_g_sq_neg_I} = {total_preimages}")
    print(f"The number of involutions is the number of preimages divided by the size of the center ({center_size}).")
    print(f"Number of involutions in PSL(4,3) = {total_preimages} / {center_size} = {inv_psl_4_3}")
    print("-" * 30)

    # --- Conclusion ---
    print(f"Result for PSL(3,9): {inv_psl_3_9} involutions.")
    print(f"Result for PSL(4,3): {inv_psl_4_3} involutions.")
    
    if inv_psl_3_9 == inv_psl_4_3:
        print("\nConclusion: PSL(3,9) and PSL(4,3) have an equal number of involutions.")
    else:
        print("\nConclusion: The groups do not have an equal number of involutions.")

calculate_involutions()