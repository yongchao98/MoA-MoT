import math

def main():
    """
    This script solves the three questions about root systems of d-neighbors of Z^n.
    It provides a step-by-step explanation for each answer based on the properties
    of the 'glue vector' w that defines the sublattice M = Z^n intersect N.
    """

    # --- Question (a) ---
    # Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?
    print("--- Analysis for Question (a) ---")
    n_a = 12
    d_a = 3
    w_a = [1] * n_a
    
    w_dot_w_a = sum(x*x for x in w_a)
    # For gcd(w_1,...,w_n, d)=1, we can check gcd of unique elements of w and d.
    g_a = math.gcd(w_a[0], d_a) 
    
    is_valid_neighbor = (w_dot_w_a % d_a == 0) and (g_a == 1)
    ans_a = "Yes" if is_valid_neighbor else "No"

    print(f"For n={n_a}, we test if R_2(M) can be A_11.")
    print(f"Let's choose d={d_a} and the glue vector w = {w_a}.")
    print(f"First, we check the neighbor conditions:")
    print(f"1. w.w must be a multiple of d. Here, w.w = {w_dot_w_a}. The equation is {w_dot_w_a} mod {d_a} = {w_dot_w_a % d_a}. This holds.")
    print(f"2. gcd(w_1,...,w_n, d) must be 1. Here, gcd(1, 3) = {g_a}. This holds.")
    print(f"Since the conditions are met, a valid neighbor lattice exists.")
    print(f"Next, we check the roots in R_2(M) = {{v in Z^12 | v.v=2, v.w = 0 mod {d_a}}}.")
    print(f"For a root v = e_i - e_j (an A_11-type root), v.w = w_i - w_j = 1 - 1 = 0. 0 mod 3 = 0. So, it's in R_2(M).")
    print(f"For a root v = e_i + e_j (a D_12-type root not in A_11), v.w = w_i + w_j = 1 + 1 = 2. 2 mod 3 != 0. So, it's not in R_2(M).")
    print(f"This construction yields R_2(M) = A_11 exactly.")


    # --- Question (b) ---
    # Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component?
    print("\n--- Analysis for Question (b) ---")
    n_b = 15
    d_b = 4
    # To create orthogonal components, we assign different values to w on different index sets.
    # Let I1 = {1..7}, I2 = {8..11}, I3 = {12..15}.
    # To isolate a D_7 component on I1, we need to make sure mixed roots are not in M.
    w_b = [0]*7 + [2]*4 + [1]*4
    
    w_dot_w_b = sum(x*x for x in w_b)
    g_b = math.gcd(math.gcd(0, 2), math.gcd(1, d_b))
    is_valid_neighbor_b = (w_dot_w_b % d_b == 0) and (g_b == 1)
    ans_b = "yes" if is_valid_neighbor_b else "no"

    print(f"For n={n_b}, we test if R_2(M) can have a D_7 component.")
    print(f"Let's choose d={d_b} and w = {w_b}.")
    print(f"First, we check the neighbor conditions:")
    print(f"1. w.w = {w_dot_w_b}. The equation is {w_dot_w_b} mod {d_b} = {w_dot_w_b % d_b}. This holds.")
    print(f"2. gcd(0, 2, 1, 4) = {g_b}. This holds.")
    print(f"Since the conditions are met, a valid neighbor lattice exists.")
    print(f"Next, we check for a D_7 component on indices I1 = {{1..7}}, where w_i=0.")
    print(f"For a root v = ±e_i ± e_j with i,j in I1, v.w = ±w_i ± w_j = ±0 ± 0 = 0. 0 mod 4 = 0. So, all D_7 roots on I1 are in R_2(M).")
    print(f"For a mixed root v = e_i + e_j with i in I1 and j in I2={{8..11}} (where w_j=2), v.w = w_i + w_j = 0 + 2 = 2. 2 mod 4 != 0. Not in R_2(M).")
    print(f"For a mixed root v = e_i + e_j with i in I1 and j in I3={{12..15}} (where w_j=1), v.w = w_i + w_j = 0 + 1 = 1. 1 mod 4 != 0. Not in R_2(M).")
    print(f"Since roots on I1 are not connected to roots outside I1, the D_7 system on I1 is an orthogonal component of R_2(M).")

    # --- Question (c) ---
    # For n = 18 and d = 5, is it possible for R_2(M) to include more than one D_n component?
    print("\n--- Analysis for Question (c) ---")
    d_c = 5
    ans_c = "no"
    print(f"For n=18, d={d_c}, we analyze the condition for having multiple D-type components.")
    print(f"Let's assume R_2(M) has a D_k component on a set of indices I.")
    print(f"This requires all roots v = ±e_i ± e_j (for i,j in I) to be in M.")
    print(f"The condition v.w = 0 mod 5 must hold. So, ±w_i ± w_j = 0 mod 5 for all i,j in I.")
    print(f"This implies two equations: (w_i + w_j) mod 5 = 0 and (w_i - w_j) mod 5 = 0.")
    print(f"Adding them gives 2*w_i = 0 mod 5. Since gcd(2, 5) = 1, this means w_i = 0 mod 5 for all i in I.")
    print(f"Now, suppose there are two D-type components: D_k1 on indices I1 and D_k2 on I2.")
    print(f"The above logic implies w_i = 0 mod 5 for all i in I1 and w_j = 0 mod 5 for all j in I2.")
    print(f"For these to be separate components, a mixed root like v = e_i + e_j (i in I1, j in I2) must NOT be in R_2(M).")
    print(f"Let's check this root's condition: v.w = w_i + w_j = 0 + 0 = 0. 0 mod 5 = 0.")
    print(f"This means the mixed root IS in R_2(M), which connects the two supposed components into a single larger one.")
    print(f"Therefore, it is not possible to have more than one D-type component.")

    print("\nFinal Answer:")
    print(f"(a) [{ans_a}]; (b) [{ans_b}]; (c) [{ans_c}].")

if __name__ == "__main__":
    main()