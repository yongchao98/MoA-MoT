def solve_lattice_questions():
    """
    Solves the three questions about root systems of d-neighbors of Z^n.
    """
    
    print("Based on the classification of root systems for d-neighbors of Z^n, the root system R_2(M) is a direct sum of components determined by d and the partition (n_0, ..., n_{d-1}) where n_k is the number of coordinates of the glue vector equal to k mod d.")
    print("The components are of type D_{n_0}, D_{n_{d/2}} (if d is even), and A_{n_k + n_{d-k} - 1} for 1 <= k < d/2.\n")

    # --- Question 1 ---
    print("--- (a) For a d-neighbor N of Z^12, can R_2(M) be of type A_11? ---\n")
    print("Answer: Yes.")
    print("Explanation:")
    print("We have n = 12. We need the root system to be a single component of type A_11.")
    print("This requires an A-type component, A_{n_k + n_{d-k} - 1}, to equal A_11. This means n_k + n_{d-k} - 1 = 11.")
    print("Thus, the sum of partitions for a pair {k, d-k} must be n_k + n_{d-k} = 12.")
    print("Since n = 12, this implies all other partitions n_j must be 0. Specifically, n_0 = 0 and (if d is even) n_{d/2} = 0, so there are no D-type components.")
    print("Let's choose d = 3. The only pair {k, d-k} with 1 <= k < d/2 is {1, 2}.")
    print("We need to satisfy n_0=0 and n_1 + n_2 = 12. A valid partition is n_0 = 0, n_1 = 1, n_2 = 11.")
    print("This corresponds to a d=3 neighbor. The resulting root system is:")
    print("A_{n_1 + n_2 - 1} = A_{1 + 11 - 1} = A_{11}.")
    answer_a = "Yes"
    print("\n")
    
    # --- Question 2 ---
    print("--- (b) Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component? ---\n")
    print("Answer: Yes.")
    print("Explanation:")
    print("We have n = 15. We need the root system to contain a D_7 component.")
    print("D-type components are either D_{n_0} or D_{n_{d/2}} (for d even).")
    print("Let's try to obtain a D_{n_0} component. This requires setting n_0 = 7.")
    print("The remaining n - n_0 = 15 - 7 = 8 coordinates must be distributed among n_k for k > 0.")
    print("Let's choose d = 2. The remaining 8 coordinates must have the value 1 (mod 2), so n_1 = 8.")
    print("The partition is n_0 = 7, n_1 = 8. The sum is 7 + 8 = 15 = n.")
    print("For d=2 (even), d/2 = 1. The components are D_{n_0} and D_{n_{d/2}} = D_{n_1}.")
    print(f"This configuration for d=2 gives the root system D_{n_0} \u2295 D_{n_1} = D_7 \u2295 D_8.")
    print("This root system contains a D_7 component.")
    answer_b = "yes"
    print("\n")

    # --- Question 3 ---
    print("--- (c) For n=18 and d=5, is it possible for R_2(M) to include more than one D_n component? ---\n")
    print("Answer: No.")
    print("Explanation:")
    print("We have n = 18 and d = 5.")
    print("The possible D-type components of the root system R_2(M) are D_{n_0} and D_{n_{d/2}}.")
    print("The term D_{n_{d/2}} only exists if d is even.")
    print("Since d = 5 is odd, the only possible D-type component is D_{n_0}.")
    print("Therefore, for d=5, the root system can contain at most one D-type component (specifically, D_{n_0} if n_0 > 1).")
    print("It is impossible for R_2(M) to include more than one D component.")
    answer_c = "no"
    print("\n")

    # --- Final Answer Summary ---
    final_answer = f"({', '.join(['a', 'b', 'c'])}) [{answer_a}; {answer_b}; {answer_c}]."
    print("Final Answer Summary:")
    print(final_answer)

solve_lattice_questions()