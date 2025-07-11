def solve():
    """
    This function encapsulates the reasoning explained above to determine the truthfulness of each statement.
    Based on the analysis, only statement 3 is universally true.

    The reasoning is as follows:
    1. A general inequality A = O(nd^2 + sqrt(ndB)) is derived using the Cauchy-Schwarz inequality on local subgraph counts (4-cycles containing a given edge).
    2. Each statement is tested against this general inequality. A statement A = O(B/d^x + nd^y) is true for all d-regular graphs only if a derived condition on the exponents, d^(y-1-x) >= C (for a universal constant C), holds for all d >= 1.
    3. We check this exponent condition for each of the six statements.
       - Statement 1: y=2, x=0.5 -> d^0.5 >= C. Fails for small d. (N)
       - Statement 2: y=2, x=0.99 -> d^0.01 >= C. Fails for small d. (N)
       - Statement 3: y=2, x=1 -> d^0 >= C (i.e., 1 >= C). Holds for C<=1. (Y)
       - Statement 4: y=2, x=1.5 -> d^-0.5 >= C. Fails for large d. (N)
       - Statement 5: y=2.6, x=1.3 -> d^0.3 >= C. Fails for small d. (N)
       - Statement 6: y=1.5, x=1 -> d^-0.5 >= C. Fails for large d. (N)
    
    This leads to the result NNYNNN.
    """
    
    # Based on the derivation in the text block.
    results = ['N', 'N', 'Y', 'N', 'N', 'N']
    
    print("".join(results))

solve()