import numpy as np

def solve_lattice_questions():
    """
    Solves three questions about root systems of d-neighbors of integer lattices.
    """

    print("--- Question 1: Is it true that for a d-neighbor N of Z^12, R2(M) can be of type A_11? ---")
    
    # Let's test the construction.
    # We choose w = (1, 1, ..., 1) and d=3.
    # M = {v in Z^12 | v.w = 0 (mod 3)}
    # The roots of Z^12 are v = +/- e_i +/- e_j.
    # The root system A_11 can be represented by roots {e_i - e_j | 1 <= i,j <= 12, i != j}.
    
    n_a = 12
    d_a = 3
    w_a = np.ones(n_a, dtype=int)
    
    # We need to check which vectors v with v.v=2 lie in M.
    # v.w = (+/- 1 +/- 1) % d == 0
    # Let's test a root of type e_i - e_j (from A_11)
    # Pick i=1, j=2. v = e_1 - e_2. v.w = w_1 - w_2
    val1 = w_a[0] - w_a[1]
    print(f"For a root v = e_1 - e_2 from A_11, v.w = w_1 - w_2.")
    print(f"With w = ({w_a[0]}, {w_a[1]}, ...), this is {w_a[0]} - {w_a[1]} = {val1}.")
    print(f"Checking the condition: {val1} mod {d_a} = {val1 % d_a}. This is 0, so these roots are in R2(M).")

    # Let's test a root of type e_i + e_j
    # Pick i=1, j=2. v = e_1 + e_2. v.w = w_1 + w_2
    val2 = w_a[0] + w_a[1]
    print(f"\nFor a root v = e_1 + e_2, v.w = w_1 + w_2.")
    print(f"With w = ({w_a[0]}, {w_a[1]}, ...), this is {w_a[0]} + {w_a[1]} = {val2}.")
    print(f"Checking the condition: {val2} mod {d_a} = {val2 % d_a}. This is not 0, so these roots are NOT in R2(M).")

    print("\nConclusion for (a): By choosing d=3 and w=(1,...,1), R2(M) is precisely the set of roots of type A_11.")
    answer_a = "Yes"
    print(f"Answer (a): {answer_a}")


    print("\n\n--- Question 2: Can R2(M) of a d-neighbor of Z^15 contain a D_7 component? ---")
    # We want to form a D_7 component on indices I = {1, ..., 7}.
    # This requires that for i,j in I, +/-w_i +/- w_j = 0 (mod d).
    # This holds if w_i are all congruent to c mod d, and 2c = 0 (mod d).
    # To make it a component, for i in I, j not in I, +/-w_i +/- w_j != 0 (mod d).
    # Let's try d=2. The condition 2c = 0 (mod 2) is always true.
    # Let c=1. So w_i = 1 (mod 2) for i in {1..7}.
    # Let w_j = 0 (mod 2) for j in {8..15}.
    
    n_b = 15
    d_b = 2
    w_b = np.array([1]*7 + [0]*8)
    print(f"Let d={d_b} and w={w_b}")
    print("We test if this creates a D_7 component on the first 7 coordinates.")

    # Test an internal root for the D_7 part (e.g., e_1 + e_2)
    w_i, w_j = w_b[0], w_b[1]
    val_in = w_i + w_j
    print(f"\nFor a root v = e_1 + e_2 inside the D_7 part (i,j <= 7):")
    print(f"v.w = w_1 + w_2 = {w_i} + {w_j} = {val_in}.")
    print(f"Checking condition: {val_in} mod {d_b} = {val_in % d_b}. This is 0, so these roots are in R2(M).")

    # Test a connecting root (e.g., e_1 + e_8)
    w_i, w_j = w_b[0], w_b[7]
    val_out = w_i + w_j
    print(f"\nFor a connecting root v = e_1 + e_8 (i <= 7, j > 7):")
    print(f"v.w = w_1 + w_8 = {w_i} + {w_j} = {val_out}.")
    print(f"Checking condition: {val_out} mod {d_b} = {val_out % d_b}. This is NOT 0, so these roots are excluded.")
    
    print("\nConclusion for (b): This construction shows R2(M) = D_7 + D_8, which contains a D_7 component.")
    answer_b = "Yes"
    print(f"Answer (b): {answer_b.lower()}")
    

    print("\n\n--- Question 3: For n=18, d=5, is it possible for R2(M) to include more than one D_n component? ---")
    d_c = 5
    print(f"Here, n=18 and d={d_c}.")
    print("A D_k component requires w_i = c (mod d) for indices i in that component, where 2*c = 0 (mod d).")
    
    # Find all possible values of c
    possible_c = []
    print(f"We search for integers c in [0, ..., {d_c-1}] satisfying 2*c mod {d_c} == 0:")
    for c in range(d_c):
        print(f"Testing c = {c}: 2 * {c} = {2*c}. {2*c} mod {d_c} = {(2*c)%d_c}")
        if (2*c) % d_c == 0:
            possible_c.append(c)
    
    print(f"\nThe only value for c that works modulo {d_c} is c = {possible_c[0]}.")
    
    print("\nNow, suppose we have two D_k components, one with c1 and one with c2.")
    print("To be separate components, they must satisfy c1 +/- c2 != 0 (mod d).")
    c1 = possible_c[0]
    c2 = possible_c[0]
    
    print(f"Since there is only one possible value for c ({c1}), any two D_k components must use c1 = {c1} and c2 = {c2}.")
    print(f"Let's check the separation condition for c1={c1}, c2={c2}, d={d_c}:")
    val_sum = c1 + c2
    print(f"Checking c1 + c2: {c1} + {c2} = {val_sum}. {val_sum} mod {d_c} = {val_sum % d_c}.")
    if (val_sum % d_c) == 0:
        print("This is 0, which violates the separation condition (c1+c2 mod d != 0).")
    
    print("\nConclusion for (c): If two sets of indices have w_i = 0 (mod 5), the roots connecting them are also in M.")
    print("This means they are not separate components, but part of one larger D_k component.")
    print("Therefore, it's not possible to have more than one D_k component.")
    answer_c = "No"
    print(f"Answer (c): {answer_c.lower()}")

    return answer_a, answer_b, answer_c

if __name__ == '__main__':
    a, b, c = solve_lattice_questions()
    final_answer = f"(a) [{a}]; (b) [{b.lower()}]; (c) [{c.lower()}]."
    # The final answer format for the platform
    # print(f"\n\nFinal Answer String:\n{final_answer}")
    # print(f"<<<{final_answer}>>>")