import math

def solve_similarity_questions():
    """
    Solves the three questions about the similarity of diagonal matrices.
    """

    # Part (a): Is it true that A and B are similar if and only if the
    # multiplicities of each eigenvalue are identical?
    #
    # Explanation: Two matrices are similar if and only if they represent the same linear
    # operator under different bases, which implies they have the same invariants, including
    # the characteristic polynomial and minimal polynomial. For diagonalizable matrices (like
    # diagonal matrices), this is equivalent to having the same set of eigenvalues with
    # the same multiplicities. Two diagonal matrices with the same multiset of diagonal
    # entries are similar via a permutation matrix. Thus, the statement is true.
    answer_a = "Yes"

    # Part (b): For n = 3 and eigenvalues alpha, beta, gamma in F, how many
    # similarity classes exist if alpha, beta, and gamma are distinct?
    #
    # Explanation: A similarity class of a diagonal matrix is determined by the multiset of its
    # eigenvalues. We need to find the number of multisets of size 3 (for an n=3 matrix) that can
    # be formed from a set of 3 distinct elements {alpha, beta, gamma}. This is a problem
    # of combinations with repetition. The formula is C(k+n-1, k), where k is the number of items
    # to choose (the size of the multiset, k=3) and n is the number of distinct choices (n=3).
    # Number of classes = C(3 + 3 - 1, 3) = C(5, 3).
    k = 3  # dimension of the matrix (size of the multiset)
    n = 3  # number of distinct eigenvalues to choose from
    
    # Calculate C(n+k-1, k)
    # The numbers in the final equation are n=3 and k=3.
    # The calculation is C(3+3-1, 3) = C(5, 3)
    num_classes_b = math.comb(n + k - 1, k)
    answer_b = num_classes_b

    # Part (c): Does the number of similarity classes for diagonal matrices in M_n(F)
    # grow exponentially with n for fixed |F| = q?
    #
    # Explanation: The number of similarity classes is the number of multisets of size n
    # that can be formed from q elements. This is C(n + q - 1, n).
    # C(n + q - 1, n) = (n+q-1)! / (n! * (q-1)!) = [(n+q-1)(n+q-2)...(n+1)] / (q-1)!
    # For a fixed q, this is a polynomial in n of degree q-1.
    # Polynomial growth is much slower than exponential growth (e.g., a^n for a > 1).
    # Therefore, the growth is not exponential.
    answer_c = "No"

    # Print the final answer in the required format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_similarity_questions()
<<< (a) Yes; (b) 10; (c) No >>>