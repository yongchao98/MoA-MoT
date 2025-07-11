import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def get_root_lattice_det(rtype, rank):
    """Returns the determinant of a given irreducible root lattice."""
    if rtype == 'A':
        return rank + 1
    elif rtype == 'D':
        if rank < 2: return None  # D_n is defined for n >= 2
        return 4
    elif rtype == 'E':
        if rank == 6: return 3
        if rank == 7: return 2
        if rank == 8: return 1
    return None

def solve_and_print():
    """Solves the three questions and prints the results."""

    # --- Question 1 ---
    print("--- Analysis for Question 1 ---")
    # For a d-neighbor of Z^12, can R_2(M) be A_11?
    det_A11 = get_root_lattice_det('A', 11)
    is_sq_A11 = is_perfect_square(det_A11)
    print("For R_2(M) to be A_11, the lattice Q(A_11) must be a sublattice of Z^12.")
    print("This requires the determinant of Q(A_11) to be a perfect square.")
    print(f"The determinant of Q(A_11) is 11 + 1 = {det_A11}.")
    print(f"Is {det_A11} a perfect square? {is_sq_A11}.")
    print("Since 12 is not a perfect square, this is impossible.")
    ans_a = "No"

    # --- Question 2 ---
    print("\n--- Analysis for Question 2 ---")
    # For a d-neighbor of Z^15, can R_2(M) contain a D_7 component?
    # This means the full root system R could be D_7 ⊕ R', where rank(D_7)+rank(R') <= 15.
    det_D7 = get_root_lattice_det('D', 7)
    print("Let the root system R contain a D_7 component, so R = D_7 ⊕ R'.")
    print("The determinant of Q(R) must be a perfect square.")
    print(f"det(Q(R)) = det(Q(D_7)) * det(Q(R')) = {det_D7} * det(Q(R')).")
    print(f"Since det(Q(D_7)) = 4 is a perfect square, det(Q(R)) is a perfect square if and only if det(Q(R')) is.")
    print("We need to find if such an R' exists with rank <= 15 - 7 = 8.")
    print("Let's test R' = A_3 (rank 3). det(Q(A_3)) = 3 + 1 = 4, which is a perfect square.")
    print("So, a root system like D_7 ⊕ A_3 (total rank 7+3=10 <= 15) is possible.")
    ans_b = "yes"
    
    # --- Question 3 ---
    print("\n--- Analysis for Question 3 ---")
    # For n=18, d=5, can R_2(M) include more than one D_n component?
    det_D_part = get_root_lattice_det('D', 2) * get_root_lattice_det('D', 2) # Example with D_2
    print("Let the root system contain D_k1 ⊕ D_k2. The determinant of this part is 4 * 4 = 16.")
    print(f"det(Q(D_k1 ⊕ D_k2)) = {det_D_part}, which is a perfect square.")
    print("By the same logic as in Q2, the full root system R = D_k1 ⊕ D_k2 ⊕ R' is possible if Q(R') has a square determinant.")
    print("We need k1 + k2 + rank(R') <= 18.")
    print("A simple combination like D_2 ⊕ D_2 (total rank 4 <= 18) satisfies the determinant rule.")
    print("A construction for d=5 is also possible. Therefore, this is possible.")
    ans_c = "yes"

    # --- Final Answer ---
    final_answer_string = f"(a) [{ans_a}]; (b) [{ans_b}]; (c) [{ans_c}]"
    print(f"\nFinal formatted answer:")
    print(f"<<<{final_answer_string}>>>")

if __name__ == '__main__':
    solve_and_print()