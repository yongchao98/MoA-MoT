import math

def product(iterable):
    """Computes the product of elements in an iterable."""
    p = 1
    for i in iterable:
        p *= i
    return p

def gl_order(n, q):
    """Calculates the order of the general linear group GL(n,q)."""
    if n == 0:
        return 1
    order = q**(n * (n - 1) // 2) * product(q**i - 1 for i in range(1, n + 1))
    return order

def sl_order(n, q):
    """Calculates the order of the special linear group SL(n,q)."""
    if n == 0:
        return 1
    if q == 0:
        return 0
    return gl_order(n, q) // (q - 1)

def calculate_psl39_involutions():
    """
    Calculates the number of involutions in PSL(3,9).
    For q > 2 odd and n=3, PSL(3,q) = SL(3,q) if gcd(n, q-1) = 1.
    Here n=3, q=9, gcd(3, 8) = 1.
    Involutions t in SL(3,q) are conjugate to diag(1, -1, -1).
    The number of such involutions is the size of this conjugacy class.
    Size = |GL(3,9)| / |GL(1,9) x GL(2,9)|
    """
    n, q = 3, 9
    
    # Calculate order of GL(n,q) groups needed for the formula
    order_gl39 = gl_order(n, q)
    order_gl29 = gl_order(n - 1, q)
    order_gl19 = gl_order(1, q)
    
    # The centralizer in GL(3,9) is GL(1,9) x GL(2,9)
    centralizer_order = order_gl19 * order_gl29
    
    # Number of involutions is the size of the conjugacy class in GL(3,9)
    # which are all in SL(3,9)
    num_involutions = order_gl39 // centralizer_order
    return num_involutions

def calculate_psl43_involutions():
    """
    Calculates the number of involutions in PSL(4,3).
    PSL(4,3) = SL(4,3) / {I, -I}.
    Involutions in PSL(4,3) come from two sources in SL(4,3):
    1. Elements t such that t^2 = I (real involutions in SL(4,3)).
    2. Elements g such that g^2 = -I.
    """
    n, q = 4, 3
    order_sl43 = sl_order(n, q)

    # Source 1: t^2 = I, t != I, -I.
    # These are conjugate to diag(1,1,-1,-1).
    # Centralizer in GL(4,3) is GL(2,3) x GL(2,3).
    # Centralizer in SL(4,3) is C_GL / (q-1).
    order_gl23 = gl_order(2, q)
    c1_sl_order = (order_gl23 * order_gl23) // (q - 1)
    num_t_type = order_sl43 // c1_sl_order
    # Each pair {t, -t} gives one involution in PSL(4,3).
    involutions_from_t = num_t_type // 2

    # Source 2: g^2 = -I.
    # The centralizer in GL(4,3) is isomorphic to GL(2, q^2) = GL(2,9).
    # Centralizer in SL(4,3) is C_GL / (q-1).
    order_gl29 = gl_order(2, q**2)
    c2_sl_order = order_gl29 // (q - 1)
    num_g_type = order_sl43 // c2_sl_order
    # Each pair {g, -g} gives one involution in PSL(4,3).
    involutions_from_g = num_g_type // 2
    
    total_involutions = involutions_from_t + involutions_from_g
    return total_involutions

if __name__ == "__main__":
    num_psl39 = calculate_psl39_involutions()
    num_psl43 = calculate_psl43_involutions()

    print("Calculating the number of involutions for the groups in option B:")
    print(f"Number of involutions in PSL(3,9): {num_psl39}")
    print(f"Number of involutions in PSL(4,3): {num_psl43}")

    if num_psl39 == num_psl43:
        print("\nThe groups PSL(3,9) and PSL(4,3) have an equal number of involutions.")
        print(f"{num_psl39} = {num_psl43}")
    else:
        print("\nThe number of involutions is not equal.")
