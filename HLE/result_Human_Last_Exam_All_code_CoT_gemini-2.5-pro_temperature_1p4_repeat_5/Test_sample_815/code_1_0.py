def psl_3_9_involutions():
    """
    Calculates the number of involutions in PSL(3,9).
    """
    q = 9
    num_involutions = q**2 * (q**2 + q + 1)
    
    print("--- Calculation for PSL(3,9) ---")
    print("For PSL(3,9), the number of involutions is given by the formula: q^2 * (q^2 + q + 1)")
    print(f"Substituting q = {q}:")
    print(f"  {q}^2 * ({q}^2 + {q} + 1)")
    print(f"= {q**2} * ({q**2} + {q} + {1})")
    print(f"= {q**2} * {q**2 + q + 1}")
    print(f"= {num_involutions}")
    return num_involutions

def psl_4_3_involutions():
    """
    Calculates the number of involutions in PSL(4,3).
    """
    print("\n--- Calculation for PSL(4,3) ---")
    print("Involutions in PSL(4,3) correspond to elements g in SL(4,3) where g^2 is I or -I.")

    # Group orders
    q = 3
    # |GL(4,3)| = (3^4-1)(3^4-3)(3^4-9)(3^4-27)
    gl_4_3 = (q**4 - 1) * (q**4 - q) * (q**4 - q**2) * (q**4 - q**3)
    sl_4_3 = gl_4_3 // (q - 1)

    # Case 1: g^2 = I
    # Centralizer size calculation
    sl_2_3 = ((q**2 - 1) * (q**2 - q)) // (q - 1)
    c_sl_1 = sl_2_3**2 + sl_2_3**2
    N1 = sl_4_3 // c_sl_1
    print("\n1. Preimages g where g^2 = I:")
    print(f"   The number of such elements (N1) is |SL(4,3)| / |C(g_1)|")
    print(f"   |SL(4,3)| = {sl_4_3}")
    print(f"   |C(g_1)| = {c_sl_1}")
    print(f"   N1 = {sl_4_3} / {c_sl_1} = {N1}")

    # Case 2: g^2 = -I
    # Centralizer size calculation
    q2 = q**2
    sl_2_9 = ((q2**2 - 1) * (q2**2 - q2)) // (q2 - 1)
    c_sl_3 = 4 * sl_2_9
    N3 = sl_4_3 // c_sl_3
    print("\n2. Preimages g where g^2 = -I:")
    print(f"   The number of such elements (N3) is |SL(4,3)| / |C(g_3)|")
    print(f"   |SL(4,3)| = {sl_4_3}")
    print(f"   |C(g_3)| = {c_sl_3}")
    print(f"   N3 = {sl_4_3} / {c_sl_3} = {N3}")

    # Total involutions in PSL(4,3)
    total_preimages = N1 + N3
    num_involutions = total_preimages // 2
    print("\nTotal number of involutions in PSL(4,3) is (N1 + N3) / 2:")
    print(f"= ({N1} + {N3}) / 2")
    print(f"= {total_preimages} / 2")
    print(f"= {num_involutions}")

    return num_involutions

if __name__ == '__main__':
    t_psl_3_9 = psl_3_9_involutions()
    t_psl_4_3 = psl_4_3_involutions()

    print("\n--- Comparison ---")
    print(f"Number of involutions in PSL(3,9): {t_psl_3_9}")
    print(f"Number of involutions in PSL(4,3): {t_psl_4_3}")

    if t_psl_3_9 == t_psl_4_3:
        print("\nThe number of involutions is equal for PSL(3,9) and PSL(4,3).")
    else:
        print("\nThe number of involutions is NOT equal.")
