def solve_lattice_questions():
    """
    This script explains the solutions to the three questions about lattice root systems
    by demonstrating constructive proofs for each case.
    """

    # --- Part 1 ---
    print("--- Question 1: Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11? ---")
    print("Plan: We will construct a sublattice M of Z^12 of index d and show its root system R_2(M) is A_11.")
    print("The existence of the corresponding d-neighbor N is guaranteed by lattice theory.\n")

    n_a = 12
    # Define a vector w in Z^12. The simplest non-trivial choice is the all-ones vector.
    w_a = [1] * n_a
    # The index 'd' is determined by the norm of w.
    d_a = sum(x*x for x in w_a)
    
    print(f"1. Let L = Z^{n_a}. We choose a primitive vector w = {tuple(w_a)} in L.")
    print(f"2. The squared norm of w is w.w = {d_a}. We will work with d = {d_a}.")
    print(f"3. We define a sublattice M = {{x in L | w.x is divisible by d}}.")
    print(f"   The index of M in L is [L:M] = {d_a}.")
    print("4. We check which vectors of norm 2 (roots) of L are in M. These roots are of the form e_i-e_j or e_i+e_j.")

    # For a root x = e_i - e_j, its dot product with w is w_i - w_j.
    dot_prod_type1 = w_a[0] - w_a[1]
    # For a root x = e_i + e_j, its dot product with w is w_i + w_j.
    dot_prod_type2 = w_a[0] + w_a[1]

    print(f"   - For a root x = e_i - e_j, w.x = w_i - w_j = 1 - 1 = {dot_prod_type1}.")
    print(f"     Since {dot_prod_type1} is divisible by {d_a}, all roots of the form e_i-e_j are in M.")
    print(f"   - For a root x = e_i + e_j, w.x = w_i + w_j = 1 + 1 = {dot_prod_type2}.")
    print(f"     Since {dot_prod_type2} is not divisible by {d_a}, these roots are not in M.")

    print("\n5. The set of roots in M, R_2(M), is therefore {+-(e_i - e_j) | 1<=i<j<=12}.")
    print("   This set of vectors forms a root system of type A_11.")
    print("Conclusion: Yes, it is possible.")
    answer_a = "Yes"

    # --- Part 2 ---
    print("\n--- Question 2: Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component? ---")
    print("Plan: We show that for d=2, the visible part M can be the D_15 lattice, whose root system R_2(M) is D_15.")
    print("We then show that the D_15 root system contains D_7 as a direct sum component.\n")
    
    n_b = 15
    d_b = 2
    print(f"1. Let L = Z^{n_b}. Consider the sublattice M consisting of all vectors with an even sum of coordinates.")
    print(f"   This is the D_{n_b} lattice. Its index in L is [L:M] = {d_b}.")
    print(f"2. A d-neighbor for d={d_b} exists giving this intersection M.")
    print("3. The roots of L are vectors x = +-(e_i +- e_j). The sum of coordinates for such a vector is -2, 0, or 2, all of which are even.")
    print("   Therefore, all roots of Z^15 are contained in M. So, R_2(M) is the full D_15 root system.")
    print("4. A root system D_k can be decomposed into orthogonal components. The root system D_15, built on basis vectors {e_1,...,e_15}}, can be written as a direct sum:")
    print("   D_15 = D_7 \u2295 D_8")
    print("   where D_7 is built on {e_1,...,e_7} and D_8 is built on {e_8,...,e_15}.")
    print("Conclusion: Since R_2(M) can be D_15, it can contain a D_7 component. So, yes.")
    answer_b = "yes"
    
    # --- Part 3 ---
    print("\n--- Question 3: For n=18, d=5, is it possible for R_2(M) to include more than one D_n component? ---")
    print("Plan: We construct a sublattice M for n=18 and d=5 and show its root system contains a direct sum of D_k systems.\n")
    
    n_c = 18
    d_c = 5
    # Choose a primitive vector w in Z^18 with w.w = 5.
    w_c = [1, 1, 1, 1, 1] + [0] * (n_c - d_c)
    
    print(f"1. Let L = Z^{n_c} and d={d_c}. We choose a primitive vector w in L with w.w={d_c}.")
    print(f"   An example is w = {tuple(w_c)}.")
    print(f"2. We construct M = {{x in L | w.x is divisible by {d_c}}}.")
    print(f"3. We check which roots x = +-(e_i +- e_j) of L are in M.")
    
    print("   Consider roots built on basis vectors {e_6, ..., e_18}, where w has zero components.")
    print("   For any such root x, its dot product w.x will be 0, which is divisible by 5.")
    print("   Therefore, all roots of the D_13 system built on {e_6, ..., e_18} are in R_2(M).")
    
    num_zero_coords = n_c - d_c
    k = 2
    l = num_zero_coords - k
    print(f"4. The root system D_{num_zero_coords} itself contains orthogonal D-type components.")
    print(f"   For example, we can decompose it as a direct sum: D_{num_zero_coords} = D_{k} \u2295 D_{l}, where {num_zero_coords} = {k} + {l}.")
    print(f"   Let's take k=2, l=11. The roots of D_2 (on {e_6,e_7}) and D_11 (on {e_8,...,e_18}) are all in M.")
    print(f"Conclusion: R_2(M) contains D_2 \u2295 D_11, which has more than one D_n component. So, yes.")
    answer_c = "yes"
    
    # --- Final Answer ---
    print("\n--- Summary of Answers ---")
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}.")

if __name__ == '__main__':
    solve_lattice_questions()
