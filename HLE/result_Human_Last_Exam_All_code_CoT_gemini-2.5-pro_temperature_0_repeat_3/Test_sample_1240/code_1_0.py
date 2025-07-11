def solve_lattice_questions():
    """
    This script analyzes three questions about the root systems of d-neighbors of Z^n
    and prints a detailed explanation and the final answers.
    """

    # --- Part (a) ---
    print("Part (a): Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?")
    print("Plan: We will attempt to construct such a neighbor lattice.")
    print("1. The root system of Z^12, R_2(Z^12), is of type D_12. Its roots are vectors v = +/- e_i +/- e_j for 1 <= i < j <= 12.")
    print("2. The visible root system R_2(M) is the subset of these D_12 roots that lie in the sublattice M = Z^12 \cap N.")
    print("3. The root system A_11 can be represented by the set of vectors {+/-(e_i - e_j) | 1 <= i < j <= 12}.")
    print("4. We need to find a sublattice M of Z^12 such that its roots of norm 2 are exactly this set.")
    print("5. Let's define M using a congruence relation based on a vector w and a modulus d. Let w = (1, 1, ..., 1) and d = 3.")
    print("   M = {x in Z^12 | w . x = sum(x_i) === 0 (mod 3)}. This M is a sublattice of index 3.")
    print("6. We check which roots of D_12 are in M by checking if their coordinates sum to a multiple of 3:")
    print("   - For a root v = e_i - e_j, the sum of coordinates is 1 + (-1) = 0. Since 0 % 3 == 0, this root is in M.")
    print("   - For a root v = e_i + e_j, the sum of coordinates is 1 + 1 = 2. Since 2 % 3 != 0, this root is NOT in M.")
    print("   - The negative roots -(e_i - e_j) and -(e_i + e_j) follow the same logic, with sums 0 and -2 respectively.")
    print("7. The set of roots in M is {+/-(e_i - e_j)}, which is precisely the root system A_11.")
    answer_a = "Yes"
    print(f"Conclusion for (a): {answer_a}, it is possible.\n")

    # --- Part (b) ---
    print("Part (b): Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component?")
    print("Plan: We will construct a lattice M whose root system contains a D_7 component.")
    print("1. A D_7 root system is formed by the vectors {+/- e_i +/- e_j | 1 <= i < j <= 7}.")
    print("2. We need to find a sublattice M of Z^15 that contains all these vectors.")
    print("3. Let's define M with d = 2 and a vector w that has 1s for the first 7 components and 0s for the rest.")
    print("   w = (1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0).")
    print("   M = {x in Z^15 | w . x = x_1 + ... + x_7 === 0 (mod 2)}. This is a sublattice of index 2.")
    print("4. Let's check if the D_7 roots v = +/- e_i +/- e_j (with 1 <= i < j <= 7) are in M.")
    print("   The dot product w . v is w_i +/- w_j = 1 +/- 1. This results in either 0 or 2.")
    print("   - Equation for one case: 1 + 1 = 2. Since 2 % 2 == 0, the root is in M.")
    print("   - Equation for another case: 1 - 1 = 0. Since 0 % 2 == 0, the root is in M.")
    print("5. All vectors of the D_7 system are in M. Therefore, R_2(M) contains a D_7 component.")
    answer_b = "yes"
    print(f"Conclusion for (b): {answer_b}.\n")

    # --- Part (c) ---
    print("Part (c): For n = 18 and d = 5, is it possible for R_2(M) to include more than one D_n component?")
    print("Plan: We will analyze the conditions required for a D_k component to exist in R_2(M).")
    print("1. Let M = {x in Z^18 | w . x === 0 (mod 5)} for some primitive vector w.")
    print("2. For R_2(M) to contain a D_k component on a set of indices I_k, all roots {+/- e_i +/- e_j | i, j in I_k} must be in M.")
    print("3. This requires w . (e_i + e_j) = w_i + w_j === 0 (mod 5) and w . (e_i - e_j) = w_i - w_j === 0 (mod 5).")
    print("4. Adding these two congruences gives the equation: (w_i + w_j) + (w_i - w_j) = 2*w_i === 0 (mod 5).")
    print("5. Since 5 is prime and does not divide 2, this implies w_i === 0 (mod 5).")
    print("6. So, to have a D_k component on indices I_k, we must have w_i === 0 (mod 5) for all i in I_k.")
    print("7. Now, assume R_2(M) contains two disjoint components, D_k on indices I_k and D_l on indices I_l.")
    print("   This requires w_i === 0 (mod 5) for all i in I_k and w_j === 0 (mod 5) for all j in I_l.")
    print("8. For the components to be separate, 'cross roots' like v = e_i + e_j (with i in I_k, j in I_l) must NOT be in M.")
    print("9. Let's check the condition for this cross root: w . v = w_i + w_j.")
    print("   The equation is: w_i + w_j === 0 + 0 = 0 (mod 5).")
    print("10. This means the cross root IS in M, which contradicts the components being separate. The roots on I_k and I_l connect to form a single D_{k+l} component.")
    answer_c = "no"
    print(f"Conclusion for (c): {answer_c}, it is not possible.\n")

    # --- Final Answer ---
    final_answer_string = f"({answer_a}); ({answer_b}); ({answer_c})"
    print("Final Answer:")
    print(final_answer_string)
    print(f"<<<{final_answer_string}>>>")

solve_lattice_questions()