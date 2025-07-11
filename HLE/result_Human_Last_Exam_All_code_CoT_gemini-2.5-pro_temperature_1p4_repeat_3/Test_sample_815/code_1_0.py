import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def order_psl(n, q):
    """Computes the order of PSL(n, q)."""
    d = gcd(n, q - 1)
    order = q**(n * (n - 1) // 2)
    for i in range(2, n + 1):
        order *= (q**i - 1)
    return order // d

def order_psu(n, q):
    """Computes the order of PSU(n, q) = |SU(n,q)|/d."""
    d = gcd(n, q + 1)
    order = q**(n * (n - 1) // 2)
    for i in range(2, n + 1):
        order *= (q**i - (-1)**i)
    return order // d

def involutions_psl34():
    """Calculates number of involutions in PSL(3,4)."""
    n, q = 3, 4
    # For PSL(n, q) with q even, involutions are transvections. There's one class.
    # Centralizer C of a transvection in SL(3,q) has order q^3*(q-1).
    # |Z(SL(3,q))| = gcd(3, q-1) = 3.
    # |C_PSL| = |C_SL| / |Z|
    group_order = order_psl(n, q)
    c_sl_order = q**3 * (q - 1)
    z_order = gcd(n, q - 1)
    c_psl_order = c_sl_order // z_order
    num_involutions = group_order // c_psl_order
    
    print("1. For PSL(3,4):")
    print(f"   Order of PSL(3,4) is {group_order}.")
    print(f"   Involutions are transvections, forming a single conjugacy class.")
    print(f"   The centralizer of an involution has order |C| = |C_SL|/|Z| = {c_sl_order}/{z_order} = {c_psl_order}.")
    print(f"   Number of involutions = |PSL(3,4)| / |C| = {group_order} / {c_psl_order} = {num_involutions}.")
    print("-" * 20)
    return num_involutions

def involutions_psu33():
    """Calculates number of involutions in PSU(3,3)."""
    n, q = 3, 3
    # For PSU(n, q) with n odd, q odd, involutions in SU(n,q) have eigenvalues {-1, ..., -1, 1}.
    # k must be even for det=1, so k=2 for n=3. One class of involutions.
    # PSU(3,3) = SU(3,3) since gcd(3, 3+1)=1.
    # Centralizer C of an involution in SU(n,q) has order |C_GU|/(q+1)
    # where C_GU = GU(n-1,q) x GU(1,q).
    def order_gu(k, q_):
        order = q_**(k * (k-1) // 2)
        for i in range(1, k + 1):
            order *= (q_**i - (-1)**i)
        return order
        
    group_order = order_psu(n, q)
    c_gu_order = order_gu(n - 1, q) * order_gu(1, q)
    c_su_order = c_gu_order // (q + 1)
    num_involutions = group_order // c_su_order
    
    print("2. For PSU(3,3):")
    print(f"   Order of PSU(3,3) is {group_order}.")
    print(f"   Involutions correspond to matrices with eigenvalues (-1, -1, 1).")
    print(f"   The centralizer has order |C| = (|GU(2,3)|*|GU(1,3)|)/(3+1) = {c_gu_order}/4 = {c_su_order}.")
    print(f"   Number of involutions = |PSU(3,3)| / |C| = {group_order} / {c_su_order} = {num_involutions}.")
    print("-" * 20)
    return num_involutions

def involutions_psl39():
    """Calculates number of involutions in PSL(3,9)."""
    n, q = 3, 9
    # For PSL(n, q) with n odd, q odd, Z(SL(n,q)) is trivial if gcd(n,q-1)=1.
    # gcd(3, 8) = 1. Involutions in PSL(3,9) correspond to involutions in SL(3,9).
    # Eigenvalues must be {-1, -1, 1} for det=1. One class.
    # Centralizer C in SL(n,q) is C_GL / (q-1), where C_GL = GL(n-1,q) x GL(1,q).
    def order_gl(k, q_):
        order = q_**(k*(k-1)//2)
        for i in range(1, k + 1):
            order *= (q_**i - 1)
        return order

    group_order = order_psl(n, q)
    c_gl_order = order_gl(n - 1, q) * order_gl(1, q)
    c_sl_order = c_gl_order // (q - 1)
    num_involutions = group_order // c_sl_order

    print("3. For PSL(3,9):")
    print(f"   Order of PSL(3,9) is {group_order}.")
    print(f"   Involutions have preimages in SL(3,9) with eigenvalues (-1, -1, 1).")
    print(f"   The centralizer has order |C| = (|GL(2,9)|*|GL(1,9)|)/(9-1) = {c_gl_order}/8 = {c_sl_order}.")
    print(f"   Number of involutions = |PSL(3,9)| / |C| = {group_order} / {c_sl_order} = {num_involutions}.")
    print("-" * 20)
    return num_involutions

def involutions_psl43():
    """Calculates number of involutions in PSL(4,3)."""
    n, q = 4, 3
    # In PSL(4,3), g=AZ is an involution if A^2 is in Z(SL(4,3))={I, -I}.
    # Two types of involution classes.
    group_order = order_psl(n, q)
    z_order = gcd(n, q - 1)

    # Type 1: A^2 = I. Eigenvalues {-1, -1, 1, 1}.
    def order_gl(k, q_):
        order = q_**(k*(k-1)//2)
        for i in range(1, k + 1):
            order *= (q_**i - 1)
        return order
    c_gl1_order = order_gl(2, q) * order_gl(2, q)
    # |C_SL| = |SL(2,3)|^2 + |GL(2,3)-SL(2,3)|^2 = 2 * |SL(2,3)|^2
    sl23_order = order_gl(2, q) // (q-1)
    c_sl1_order = 2 * sl23_order**2
    c_psl1_order = c_sl1_order // z_order
    num_inv1 = group_order // c_psl1_order

    # Type 2: A^2 = -I. Preimage A is conjugate to diag(J,J) where J=[[0,-1],[1,0]].
    # Centralizer in GL(4,3) is GL(2, 3^2)=GL(2,9).
    c_gl2_order = order_gl(2, q**2)
    # C_SL is kernel of Norm determinant map, has index q-1. But for PGL to PSL, it is index 2.
    c_sl2_order = c_gl2_order // (q-1) # This is C_SL
    c_psl2_order = c_sl2_order // z_order # This is C_PSL
    num_inv2 = group_order // c_psl2_order

    total_involutions = num_inv1 + num_inv2
    
    print("4. For PSL(4,3):")
    print(f"   Order of PSL(4,3) is {group_order}.")
    print(f"   It has two classes of involutions.")
    print(f"   Class 1 (A^2=I): {num_inv1} involutions (centralizer order {c_psl1_order}).")
    print(f"   Class 2 (A^2=-I): {num_inv2} involutions (centralizer order {c_psl2_order}).")
    print(f"   Total number of involutions = {num_inv1} + {num_inv2} = {total_involutions}.")
    print("-" * 20)
    return total_involutions

def involutions_psu44():
    """Calculates number of involutions in PSU(4,4)."""
    n, q = 4, 4
    # For PSU(n, q) with q even, involutions are unipotent.
    # PSU(4,4)=SU(4,4). Two classes of involutions.
    # Class sizes from K. Shinoda's paper on conjugacy classes.
    group_order = order_psu(n, q)
    # Class 1 (transvections): (q^4-1)(q^3+1)/(q^2-1)
    num_inv1 = (q**4 - 1) * (q**3 + 1) // (q**2 - 1)
    # Class 2: q^2 * (q^4-1)
    num_inv2 = q**2 * (q**4 - 1)
    total_involutions = num_inv1 + num_inv2

    print("5. For PSU(4,4):")
    print(f"   Order of PSU(4,4) is {group_order}.")
    print(f"   It has two classes of involutions (unipotent types).")
    print(f"   Class 1 size = (4^4-1)(4^3+1)/(4^2-1) = {num_inv1}.")
    print(f"   Class 2 size = 4^2 * (4^4-1) = {num_inv2}.")
    print(f"   Total number of involutions = {num_inv1} + {num_inv2} = {total_involutions}.")
    print("-" * 20)
    return total_involutions

if __name__ == '__main__':
    counts = {}
    counts['PSL(3,4)'] = involutions_psl34()
    counts['PSU(3,3)'] = involutions_psu33()
    counts['PSL(3,9)'] = involutions_psl39()
    counts['PSL(4,3)'] = involutions_psl43()
    counts['PSU(4,4)'] = involutions_psu44()

    print("\nSummary of involution counts:")
    for group, count in counts.items():
        print(f"{group}: {count}")

    print("\nComparing pairs:")
    print(f"A. PSL(3,4) [{counts['PSL(3,4)']}] vs PSU(3,3) [{counts['PSU(3,3)']}]")
    print(f"B. PSL(3,9) [{counts['PSL(3,9)']}] vs PSL(4,3) [{counts['PSL(4,3)']}]")
    print(f"C. PSL(3,9) [{counts['PSL(3,9)']}] vs PSU(4,4) [{counts['PSU(4,4)']}]")
    print(f"D. PSL(3,4) [{counts['PSL(3,4)']}] vs PSL(3,9) [{counts['PSL(3,9)']}]")
    
    if counts['PSL(3,4)'] == counts['PSU(3,3)']:
        print("\nConclusion: A has an equal number of involutions.")
    elif counts['PSL(3,9)'] == counts['PSL(4,3)']:
        print("\nConclusion: B has an equal number of involutions.")
    elif counts['PSL(3,9)'] == counts['PSU(4,4)']:
        print("\nConclusion: C has an equal number of involutions.")
    elif counts['PSL(3,4)'] == counts['PSL(3,9)']:
        print("\nConclusion: D has an equal number of involutions.")
    else:
        print("\nConclusion: None of the pairs A, B, C, D have an equal number of involutions.")
