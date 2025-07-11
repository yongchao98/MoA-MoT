def solve_matrix_similarity_questions():
    """
    This function solves the three-part question about the similarity of diagonal matrices
    and prints the answer in the specified format.
    """

    # Part (a): Is it true that two diagonal matrices A and B are similar if and only if
    # the multiplicities of each eigenvalue are identical?
    #
    # Reasoning:
    # 1. If A and B are similar, they have the same characteristic polynomial,
    #    so they have the same eigenvalues with the same multiplicities.
    # 2. If diagonal matrices A and B have the same multiset of eigenvalues, their
    #    diagonal entries are a permutation of each other. This permutation can be
    #    achieved by a similarity transformation with a permutation matrix P, such that B = PAP⁻¹.
    # Thus, the statement is an "if and only if" condition.
    answer_a = "Yes"

    # Part (b): For n = 3 and distinct eigenvalues α, β, γ, how many similarity classes exist?
    #
    # Reasoning:
    # A similarity class of a diagonal matrix is determined by its multiset of eigenvalues.
    # The problem specifies that the eigenvalues for any matrix in consideration are {α, β, γ}.
    # Since this multiset is fixed, all such diagonal matrices belong to the same, single similarity class.
    answer_b = 1

    # Part (c): Does the number of similarity classes for diagonal matrices in M_n(F)
    # grow exponentially with n for fixed q (where F is a finite field F_q)?
    #
    # Reasoning:
    # The number of similarity classes equals the number of possible multisets of n eigenvalues
    # chosen from the q elements of the field F_q.
    # This is a "stars and bars" combinatorics problem, and the number of such multisets is
    # given by the binomial coefficient C(n + q - 1, n).
    # For a fixed q, this is a polynomial in n of degree q-1.
    # Polynomial growth is not exponential growth.
    answer_c = "No"

    # Print the final result in the specified format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}.")

solve_matrix_similarity_questions()