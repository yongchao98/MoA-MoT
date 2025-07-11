import numpy as np
import math

def solve_lattice_questions():
    """
    This script provides reasoning and answers for the three questions about
    the root systems of d-neighbors of integer lattices.
    """

    # --- Question 1 ---
    print("--- Question 1: A_11 component in a neighbor of Z^12? ---")
    n_1 = 12
    d_1 = 3
    v_1 = np.ones(n_1, dtype=int)
    print(f"Plan: Construct a sublattice M of Z^{n_1} using the Kneser construction with d={d_1}.")
    print(f"Let M = {{ x in Z^n | v.x = 0 mod d }} with n={n_1}, d={d_1}, and v = (1, 1, ..., 1).")
    # Check for the existence of the d-neighbor N
    v_dot_v_1 = np.dot(v_1, v_1)
    print(f"A d-neighbor exists if v.v is divisible by d. Here v.v = {v_dot_v_1}.")
    print(f"Equation: {v_dot_v_1} mod {d_1} = {v_dot_v_1 % d_1}. The condition holds.")
    # Determine the root system R2(M)
    # Check roots of type e_i - e_j
    print("The roots of Z^12 are of type e_i +/- e_j. We test which ones are in M.")
    dot_u1 = v_1[0] - v_1[1]
    print(f"For u = e_i - e_j: v.u = {v_1[0]} - {v_1[1]} = {dot_u1}. Since {dot_u1} mod {d_1} = 0, these roots are in M.")
    # Check roots of type e_i + e_j
    dot_u2 = v_1[0] + v_1[1]
    print(f"For u = e_i + e_j: v.u = {v_1[0]} + {v_1[1]} = {dot_u2}. Since {dot_u2} mod {d_1} != 0, these roots are NOT in M.")
    print("Conclusion: R2(M) is exactly the A_11 root system.")
    print("(a) Yes")

    print("\n" + "="*50 + "\n")

    # --- Question 2 ---
    print("--- Question 2: D_7 component in a neighbor of Z^15? ---")
    n_2 = 15
    k_2 = 7 # Size of D_7 component
    d_2 = 4
    a_2 = 2
    b_2 = 1
    v_2 = np.array([a_2] * k_2 + [b_2] * (n_2 - k_2), dtype=int)
    print(f"Plan: Construct a suitable M for n={n_2} with d={d_2}.")
    print(f"We use a vector v with two parts: {k_2} entries of '{a_2}' and {n_2 - k_2} entries of '{b_2}'.")
    v_dot_v_2 = np.dot(v_2, v_2)
    print(f"Check neighbor existence: v.v = {k_2}*({a_2}^2) + {n_2 - k_2}*({b_2}^2) = {v_dot_v_2}.")
    print(f"Equation: {v_dot_v_2} mod {d_2} = {v_dot_v_2 % d_2}. The condition holds.")
    print("Check which roots of D_15 are in M:")
    # D_7 component on first k_2 indices
    dot1_p, dot1_m = a_2 + a_2, a_2 - a_2
    print(f"  For D_7 part (i,j <= 7): v.(e_i+/-e_j) = {a_2}+/-{a_2} gives {dot1_p} or {dot1_m}. Both are 0 mod {d_2}, so D_7 is in R2(M).")
    # Roots on remaining indices
    dot2_p, dot2_m = b_2 + b_2, b_2 - b_2
    print(f"  For other indices (i,j > 7): v.(e_i-e_j) = {b_2}-{b_2} = {dot2_m} is 0 mod {d_2}. v.(e_i+e_j) = {b_2}+{b_2} = {dot2_p} is not.")
    # Cross terms
    dot3_p, dot3_m = a_2 + b_2, a_2 - b_2
    print(f"  For cross terms (i<=7, j>7): v.(e_i+/-e_j) = {a_2}+/-{b_2} gives {dot3_p} or {dot3_m}. Neither is 0 mod {d_2}.")
    print("Conclusion: R2(M) is of type D_7 + A_7, which contains a D_7 component.")
    print("(b) yes")

    print("\n" + "="*50 + "\n")

    # --- Question 3 ---
    print("--- Question 3: More than one D_n component for n=18, d=5? ---")
    n_3 = 18
    d_3 = 5
    print("Plan: Proof by contradiction.")
    print("Assume R2(M) has two components, D_k on index set I_k and D_l on I_l.")
    print("For M={x|v.x=0 mod 5}, having a D_k component implies 2*v_i = 0 mod 5 for i in I_k.")
    print(f"Since 2 and 5 are coprime, this implies v_i = 0 (mod {d_3}).")
    print(f"The same holds for the D_l component: v_j = 0 (mod {d_3}) for j in I_l.")
    print("Now consider a cross-root u = e_i + e_j, with i in I_k and j in I_l.")
    print(f"v.u = v_i + v_j. Since v_i = 0 (mod 5) and v_j = 0 (mod 5), v.u = 0 + 0 = 0 (mod 5).")
    print("Conclusion: All cross-roots are present, merging D_k and D_l into a single D_{k+l} component.")
    print("The root system cannot be decomposed into more than one D-type component.")
    print("(c) no")

solve_lattice_questions()

# Present the final answer in the required format
final_answer = "(a) Yes; (b) yes; (c) no"
print(f"\n<<<{final_answer}>>>")