import numpy as np

def solve_questions():
    """
    Solves the three questions by constructing witness vectors or proving non-existence.
    """

    # --- Question 1 ---
    print("--- Question 1 ---")
    print("Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?\n")
    print("Plan: We need to find a primitive vector w in Z^12 and a d=w.w such that")
    print("the visible root system R_2(M) = {v in R_2(Z^12) | v.w == 0 (mod d)} is of type A_11.")
    
    n_a = 12
    # To obtain an A_11 system, we need a set of 12 indices I_c such that for i, j in I_c,
    # w_i = w_j (mod d). This means all 12 coordinates of w must be congruent modulo d.
    # The simplest primitive vector satisfying this is w = (1, 1, ..., 1).
    w_a = np.ones(n_a, dtype=int)
    d_a = np.dot(w_a, w_a)
    
    print(f"Let's choose the primitive vector w = {tuple(w_a)} in Z^{n_a}.")
    print(f"The squared norm is d = w . w = {' + '.join(['1^2']*n_a)} = {d_a}.")
    print(f"The sublattice M is {{x in Z^{n_a} | x.w == 0 (mod {d_a})}}.")
    
    print("\nThe roots of Z^12 are v = +/- e_i +/- e_j. Let's check which ones are in M.")
    print("1. For roots v = e_i - e_j:")
    print(f"   v.w = w_i - w_j = 1 - 1 = 0. This is always divisible by d={d_a}.")
    print("   So, all vectors of the form +/-(e_i - e_j) are in M. These form the root system A_11.")
    
    print("2. For roots v = e_i + e_j:")
    print(f"   v.w = w_i + w_j = 1 + 1 = 2. Is 2 divisible by d={d_a}? No.")
    print("   So, no vectors of the form +/-(e_i + e_j) are in M.")
    
    print("\nThe visible root system R_2(M) is therefore precisely A_11.")
    answer_a = "Yes"
    print(f"Answer: {answer_a}\n")

    # --- Question 2 ---
    print("--- Question 2 ---")
    print("Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component?\n")
    print("Plan: We need to find a primitive vector w in Z^15 and a d=w.w such that")
    print("for a set of 7 indices I_c, the roots +/- e_i +/- e_j for i,j in I_c are all in M.")
    print("This requires finding a set of 7 indices I_c where w_i == c (mod d) and 2*c == 0 (mod d).")
    
    n_b = 15
    # Let's try to satisfy the condition 2*c == 0 (mod d) with c=0.
    # This requires a set I_0 of size 7, where w_i == 0 (mod d).
    # A simple choice is to set 7 coordinates of w to 0.
    # To make the vector primitive, the other 8 coordinates can be 1.
    w_b = np.zeros(n_b, dtype=int)
    num_zeros = 7
    num_ones = n_b - num_zeros
    w_b[num_zeros:] = 1
    d_b = int(np.dot(w_b, w_b))
    
    print(f"Let's construct a vector w in Z^15 with 7 zeros and 8 ones: w = {tuple(w_b)}.")
    print(f"This vector is primitive because its entries have a greatest common divisor of 1.")
    print(f"The squared norm is d = w . w = {' + '.join(['0^2']*num_zeros)} + {' + '.join(['1^2']*num_ones)} = {d_b}.")
    
    print("\nLet's check for a D_7 component. Consider the first 7 indices I = {1, ..., 7}.")
    print("For any i in I, w_i = 0. This corresponds to the partition I_c with c=0.")
    print(f"The condition for a D component is 2*c == 0 (mod d). With c=0, d={d_b}, we have 2*0 = 0, which is divisible by {d_b}. The condition is met.")
    
    print("To verify, let i,j be in {1,...,7}. For any root v = +/- e_i +/- e_j:")
    print(f"   v.w = +/- w_i +/- w_j = +/- 0 +/- 0 = 0. This is divisible by d={d_b}.")
    print("   So all vectors +/- e_i +/- e_j for i,j in {1,...,7} are in M. These form a D_7 root system.")
    
    print("\nThus, R_2(M) contains a D_7 component.")
    answer_b = "Yes"
    print(f"Answer: {answer_b}\n")
    
    # --- Question 3 ---
    print("--- Question 3 ---")
    print("For n = 18 and d = 5, is it possible for R_2(M) to include more than one D_n component?\n")
    
    n_c, d_c = 18, 5
    
    print("Plan: We will analyze the conditions for getting a D-type component and apply it to n=18, d=5.")
    print("A D_k component arises from a partition of indices I_c = {i | w_i == c (mod d)} where k=|I_c| > 1 and 2*c == 0 (mod d).")
    print(f"For d=5, the condition becomes 2*c == 0 (mod 5).")
    print("Since 5 is prime and does not divide 2, this implies c must be a multiple of 5. So, c == 0 (mod 5).")
    
    print("\nThis means that only the partition I_0 = {i | w_i == 0 (mod 5)} can generate a D-type component.")
    print("Since there is only one such partition I_0 for any given vector w, there can be at most one D-type component.")
    
    print("\nTo confirm, we can examine all primitive vectors w in Z^18 with w.w = 5.")
    print("The only ways to write 5 as a sum of integer squares (up to permutation and sign changes) are:")
    print("1) 5 = 2^2 + 1^2")
    print("2) 5 = 1^2 + 1^2 + 1^2 + 1^2 + 1^2")

    print("\nCase 1: w has coordinates (+/-2, +/-1, 0, ..., 0). Example: w = (2, 1, 0, ...).")
    print("   w.w = 2^2 + 1^2 + 0 + ... = 5. The partitions mod 5 are:")
    print("   I_0 = {i | w_i=0} (16 indices). 2*0=0 (mod 5) -> D_16 component.")
    print("   I_1 = {i | w_i=1} (1 index).   2*1!=0 (mod 5) -> No D component.")
    print("   I_2 = {i | w_i=2} (1 index).   2*2!=0 (mod 5) -> No D component.")
    print("   The root system is D_16, which has only one D-component.")
    
    print("\nCase 2: w has five coordinates of +/-1, and the rest 0. Example: w = (1,1,1,1,1,0,...).")
    print(f"   w.w = 1^2 + 1^2 + 1^2 + 1^2 + 1^2 = 5. The partitions mod 5 are:")
    print("   I_0 = {i | w_i=0} (13 indices). 2*0=0 (mod 5) -> D_13 component.")
    print("   I_1 = {i | w_i=1} (5 indices).  2*1!=0 (mod 5) -> A_4 component.")
    print("   The root system is D_13 + A_4, which has only one D-component.")

    print("\nIn all possible cases, the resulting root system has at most one D-type component.")
    answer_c = "No"
    print(f"Answer: {answer_c}\n")

    # --- Final Answer ---
    final_answer_str = f"(a) [{answer_a}]; (b) [{answer_b.lower()}]; (c) [{answer_c.lower()}]."
    print("----------------------")
    print(final_answer_str)
    print(f"<<<{final_answer_str}>>>")

solve_questions()