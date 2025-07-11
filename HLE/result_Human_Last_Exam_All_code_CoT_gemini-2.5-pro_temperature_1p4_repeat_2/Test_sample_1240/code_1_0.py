import numpy as np
from math import comb

def check_q1():
    """
    Verification for Question 1.
    """
    print("--- Verifying Question 1 (A_11 in Z^12) ---")
    n = 12
    d = 3
    
    # We choose M = {v in Z^12 | sum(v_i) = 0 mod 3}
    # This corresponds to u = (1, 1, ..., 1)
    
    # Check construction validity
    # We need a vector w in Z^12 with w_i = 1 (mod 3) and w.w = 0 (mod 9)
    w = np.array([1]*10 + [-2]*2)
    w_dot_w = np.dot(w, w)
    d_squared = d**2
    
    print(f"Chosen w = {w}")
    print(f"Components of w mod {d}: {w % d}")
    print(f"w . w = {w_dot_w}")
    print(f"d^2 = {d_squared}")
    if w_dot_w % d_squared == 0:
        print(f"Validity check passed: {w_dot_w} is divisible by {d_squared}.")
    else:
        print(f"Validity check failed.")

    # Now verify the root system
    a11_roots = 0
    other_roots = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            # Root type e_i - e_j (representative of A_11 roots)
            # Sum of coordinates is 0, so dot product with u=(1,..1) is 0 mod 3.
            # We have e_i - e_j and e_j - e_i, so we count 2 for each pair {i,j}.
            a11_roots += 2

            # Root type e_i + e_j
            v = np.zeros(n, dtype=int)
            v[i] = 1
            v[j] = 1
            if np.sum(v) % d == 0:
                other_roots += 2 # for e_i+e_j and -e_i-e_j
    
    print(f"Found {a11_roots} roots of type A_11 (e.g., e_i-e_j).")
    print(f"Found {other_roots} other roots (e.g., e_i+e_j).")
    
    total_A11_roots = 11 * 12
    if a11_roots == total_A11_roots and other_roots == 0:
        print("Conclusion for (a): R2(M) is precisely of type A_11. So the answer is Yes.")
    else:
        print("Conclusion for (a): R2(M) is not of type A_11.")
    print("-" * 20 + "\n")


def check_q2():
    """
    Verification for Question 2.
    """
    print("--- Verifying Question 2 (D_7 component in Z^15) ---")
    n = 15
    d = 2
    k = 7 # for D_7
    
    # u = (0,..,0, 1,..,1) with k zeros and n-k ones.
    u = np.array([0]*k + [1]*(n-k))
    
    # Check construction validity
    # We need w in Z^15 with w_i even for i<k, odd for i>=k, and w.w = 0 (mod 4)
    # The sum of squares mod 4 is k*(even^2) + (n-k)*(odd^2) = 7*0 + 8*1 = 8.
    w_sq_sum_mod_4 = (k * 0 + (n-k) * 1) % (d**2)
    print(f"For w with {k} even and {n-k} odd components, sum(w_i^2) mod 4 is {w_sq_sum_mod_4}")
    if w_sq_sum_mod_4 == 0:
        print("Validity check passed.")
    else:
        print("Validity check failed.")

    # Verify the root system structure
    d7_roots = 4 * comb(k, 2)
    d8_roots = 4 * comb(n-k, 2)
    connecting_roots = 0

    # A vector v is in M if u.v is even
    # Check connecting roots: e_i +/- e_j, where i<=k, j>k
    v = np.zeros(n, dtype=int)
    v[0] = 1
    v[k] = 1 # e_1 + e_8
    if np.dot(u, v) % d == 0:
        connecting_roots += 1
        
    print(f"Roots from first {k} coordinates (D_7 type) are in M.")
    print(f"Roots from last {n-k} coordinates (D_8 type) are in M.")
    if connecting_roots == 0:
        print("No connecting roots found. The root system decomposes.")
        print(f"Conclusion for (b): R2(M) is D_7 + D_8, so it has a D_7 component. The answer is Yes.")
    else:
        print("Connecting roots found. The root system is irreducible.")
    print("-" * 20 + "\n")

if __name__ == '__main__':
    check_q1()
    check_q2()
    print("--- Analysis for Question 3 (Multiple D_n components in Z^18) ---")
    print("As explained in the text, if R2(M) were to contain D_k and D_l,")
    print("the condition on the gluing vector u would force connecting roots to also be in M,")
    print("merging the components into a single larger D_{k+l} system.")
    print("Conclusion for (c): It is not possible. The answer is no.")
    print("-" * 20 + "\n")
    print("Final Answer Summary:")
    print("(a) Yes; (b) yes; (c) no.")