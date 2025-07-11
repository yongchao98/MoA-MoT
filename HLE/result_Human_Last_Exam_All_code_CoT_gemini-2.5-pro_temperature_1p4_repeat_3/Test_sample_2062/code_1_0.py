def solve_matrix_similarity_questions():
    """
    This function provides the answers to the three questions about matrix similarity.
    """
    # (a) Is it true that A and B are similar if and only if the
    # multiplicities of each eigenvalue are identical?
    answer_a = "Yes"

    # (b) For n = 3 and distinct eigenvalues, how many similarity classes exist?
    # The number of similarity classes is 1, as the multiset of eigenvalues is fixed.
    answer_b = 1

    # (c) Does the number of similarity classes grow exponentially with n?
    # No, it grows polynomially with n. The degree of the polynomial is q-1.
    answer_c = "No"

    # The final output is formatted as requested.
    # Note for part (b), the number 1 is the expression for the answer.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}.")

solve_matrix_similarity_questions()