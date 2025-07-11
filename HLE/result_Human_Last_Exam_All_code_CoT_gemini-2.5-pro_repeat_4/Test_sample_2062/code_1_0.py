def solve_matrix_similarity_questions():
    """
    Analyzes the similarity of diagonal matrices and answers the given questions.
    """

    # Part (a): Are two diagonal matrices A and B similar if and only if
    # the multiplicities of each eigenvalue are identical?
    #
    # Reasoning: Two matrices are similar if and only if they have the same
    # Jordan Normal Form. For a diagonal matrix, its eigenvalues are its
    # diagonal entries, and it is already in Jordan form (with all Jordan blocks
    # of size 1). Two diagonal matrices have the same Jordan form if and only
    # if they have the same diagonal entries up to a permutation. This is
    # precisely the condition that they have the same eigenvalues with the
    # same multiplicities. Thus, the statement is true.
    answer_a = "Yes"

    # Part (b): For n=3 and distinct eigenvalues alpha, beta, gamma,
    # how many similarity classes exist?
    #
    # Reasoning: Any diagonal matrix with these eigenvalues must have the
    # multiset of eigenvalues {alpha, beta, gamma}, with each eigenvalue having
    # multiplicity 1. By the conclusion in part (a), all such matrices
    # are similar to each other. Therefore, there is only one similarity class.
    answer_b = 1

    # Part (c): Does the number of similarity classes for diagonal matrices in M_n(F)
    # grow exponentially with n for a fixed field size q?
    #
    # Reasoning: A similarity class is defined by the multiset of n eigenvalues
    # chosen from the field F of size q. The number of ways to do this is a
    # combinations with repetition problem, given by the formula:
    # C(n + q - 1, n) = (n + q - 1)! / (n! * (q-1)!).
    # For a fixed q, this is a polynomial in n of degree q-1.
    # Polynomial growth is not exponential. Therefore, the answer is no.
    answer_c = "No"

    # The problem asks to output the number in the final equation.
    # For part (b), the number is 1.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}.")

solve_matrix_similarity_questions()