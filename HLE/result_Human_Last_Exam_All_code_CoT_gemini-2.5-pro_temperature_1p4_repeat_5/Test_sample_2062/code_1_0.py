import math

def solve_and_print():
    """
    This function encapsulates the reasoning and prints the final answer
    to the user's question about the similarity of diagonal matrices.
    """

    # Part (a): Is it true that two diagonal matrices A and B are similar if and only if
    # the multiplicities of each eigenvalue are identical?

    # Reasoning:
    # "Only if" part: If two matrices are similar, they have the same characteristic
    # polynomial, which implies they have the same eigenvalues with the same multiplicities.
    # This is true for any pair of similar matrices, including diagonal ones.
    # "If" part: For diagonal matrices A = diag(α_1, ..., α_n) and
    # B = diag(β_1, ..., β_n), if they have the same eigenvalues with the same
    # multiplicities, it means the multiset {α_1, ..., α_n} is the same as the
    # multiset {β_1, ..., β_n}. This implies that B can be obtained from A by
    # simply reordering its diagonal elements. Such a reordering can be
    # represented by conjugation with a permutation matrix P (i.e., B = P⁻¹AP).
    # Since permutation matrices are invertible, A and B are similar.
    # Therefore, the statement is true.
    answer_a = "Yes"

    # Part (b): For n = 3 and eigenvalues α, β, γ ∈ F, how many similarity classes
    # exist if α, β, and γ are distinct?

    # Reasoning:
    # As established in part (a), a similarity class for a diagonal matrix is
    # uniquely determined by the multiset of its eigenvalues.
    # The problem states the eigenvalues are α, β, and γ, and they are distinct.
    # This means any such matrix has the eigenvalue multiset {α, β, γ}, where each
    # eigenvalue has a multiplicity of 1.
    # Since there is only one such multiset of eigenvalues, all diagonal matrices
    # with these eigenvalues belong to the same, single similarity class.
    # The number in the equation that represents the number of classes is 1.
    answer_b = 1

    # Part (c): Does the number of similarity classes for diagonal matrices in M_n(F)
    # grow exponentially with n for fixed q?

    # Reasoning:
    # A similarity class is determined by the multiset of n eigenvalues, which must be
    # chosen from the field F, where |F| = q. This is equivalent to counting the number
    # of ways to choose n items from q categories with replacement, where order doesn't matter.
    # This is a classic "stars and bars" combinatorial problem. The number of such
    # multisets (and thus the number of similarity classes) is given by the formula:
    # N(n, q) = C(n + q - 1, q - 1) = (n + q - 1)! / (n! * (q - 1)!)
    # For a fixed q, this formula is a polynomial in n of degree q - 1.
    # For example, if q = 2, N = n + 1 (linear growth). If q = 3, N = (n+2)(n+1)/2 (quadratic growth).
    # Polynomial growth is significantly slower than exponential growth (e.g., c^n for c > 1).
    # Therefore, the number of similarity classes does not grow exponentially with n.
    answer_c = "No"

    # Print the final formatted answer.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_and_print()