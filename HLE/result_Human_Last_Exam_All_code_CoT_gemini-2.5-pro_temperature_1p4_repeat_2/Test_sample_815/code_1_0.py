from functools import reduce
import math

def product(iterable):
    """Helper function to compute the product of an iterable."""
    return reduce(lambda x, y: x * y, iterable, 1)

def order_GL(n, q):
    """Computes the order of the general linear group GL(n, q)."""
    # Formula: q^(n(n-1)/2) * Product_{i=1 to n} (q^i - 1)
    if n == 0:
        return 1
    term1 = q**(n * (n - 1) // 2)
    term2 = product([(q**i - 1) for i in range(1, n + 1)])
    return term1 * term2

def get_involutions_count():
    """
    Calculates and prints the number of involutions for the groups in option B.
    """
    # --- Calculation for PSL(3, 9) ---
    n1, q1 = 3, 9
    order_sl_3_9 = order_GL(n1, q1) // (q1 - 1)
    
    # Centralizer of representative involution g=diag(-1,-1,1)
    # C_GL(g) = GL(2,q) x GL(1,q)
    order_gl_2_9 = order_GL(2, q1)
    order_gl_1_9 = q1 - 1
    order_c_gl_g1 = order_gl_2_9 * order_gl_1_9
    order_c_sl_g1 = order_c_gl_g1 // (q1 - 1)
    
    num_involutions_psl_3_9 = order_sl_3_9 // order_c_sl_g1
    
    # --- Calculation for PSL(4, 3) ---
    n2, q2 = 4, 3
    order_sl_4_3 = order_GL(n2, q2) // (q2 - 1)
    d = math.gcd(n2, q2 - 1)
    order_psl_4_3 = order_sl_4_3 // d
    
    # Case 1: g^2 = I
    # Centralizer of g=diag(-1,-1,1,1) in SL(4,3)
    order_gl_2_3 = order_GL(2, q2)
    order_c_gl_g2 = order_gl_2_3 * order_gl_2_3
    order_c_sl_g2 = order_c_gl_g2 // (q2 - 1)
    # g and -g are conjugate in SL(4,3), so |C_PSL(gZ)| = |C_SL(g)|
    order_c_psl_g2 = order_c_sl_g2
    class_size_1 = order_psl_4_3 // order_c_psl_g2

    # Case 2: g^2 = -I
    # Centralizer in GL(4,3) is isomorphic to GL(2,9)
    order_gl_2_9 = order_GL(2, 9)
    # Centralizer in SL(4,3)
    order_c_sl_h2 = order_gl_2_9 // (q2 - 1)
    # h and -h are conjugate in SL(4,3), so |C_PSL(hZ)| = |C_SL(h)|
    order_c_psl_h2 = order_c_sl_h2
    class_size_2 = order_psl_4_3 // order_c_psl_h2
    
    num_involutions_psl_4_3 = class_size_1 + class_size_2
    
    print(f"Number of involutions in PSL(3,9): {num_involutions_psl_3_9}")
    print(f"Number of involutions in PSL(4,3): {class_size_1} (type g^2=I) + {class_size_2} (type g^2=-I) = {num_involutions_psl_4_3}")
    
    if num_involutions_psl_3_9 == num_involutions_psl_4_3:
        print("\nThe number of involutions is equal for PSL(3,9) and PSL(4,3).")
    else:
        print("\nThe number of involutions is NOT equal.")

get_involutions_count()
