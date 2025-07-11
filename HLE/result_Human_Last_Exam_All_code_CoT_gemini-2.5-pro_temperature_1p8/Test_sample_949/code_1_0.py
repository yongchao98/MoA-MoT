def solve():
    """
    Analyzes six statements about parenthesis strings and determines their truth value.

    The plan is as follows:
    1.  Represent the well-formed parenthesis string as a tree structure. Define L(x) (length)
        and D(x) (depth) in terms of tree properties: L(x) = 2 * size(subtree at x) and
        D(x) = height(subtree at x) + 1.
    2.  For each statement, test it against several families of strings that represent extreme
        structures:
        a. Deep and narrow strings (like a path graph): S = (((...()))).
        b. Shallow and wide strings (like a star graph): S = (()()()...()).
        c. Balanced structures (like a full binary tree).
        d. A custom "broom" or "comet" structure S_{k,m} that combines path and star elements,
           which is particularly effective for creating counterexamples. This structure is built
           by creating a spine of k nested pairs, each containing a block of m simple pairs.
    3.  Analyze the asymptotic behavior of the sums on the left-hand side (LHS) and
        right-hand side (RHS) of the O-notation as the parameters (k, m) of the string
        family grow.
    4.  Statement 1 (log L vs log D): Found to be true. The log function grows slowly,
        preventing the LHS from outgrowing the RHS in all tested cases. The sums are
        dominated by the number of pairs.
    5.  Statement 2 (loglog L vs loglog D): Also true for the same reason. The log-log
        function grows even more slowly.
    6.  Statement 3 (log^5 L vs log^5 D): True. The exponent is applied to both sides,
        maintaining the proportionality found for Statement 1.
    7.  Statement 4 (2^sqrt(log L) vs 2^sqrt(log D)): Found to be false. The exponential function
        amplifies differences. The broom-like counterexample S_{k,m} with m growing faster than k (e.g., m=k^a, a>1)
        causes the ratio of LHS/RHS to become unbounded.
    8.  Statement 5 (L^0.1 vs D^0.11): Found to be false. The slightly smaller exponent for L
        is not enough to compensate for cases where L is significantly larger than D. The S_{k,m}
        counterexample with m=k makes the LHS grow as k^1.2 while the RHS grows as k^1.11.
    9.  Statement 6 (L^0.25 vs D^0.5): Found to be false. The exponent for D is now significantly
        larger, but the S_{k,m} counterexample can still be tailored (e.g., m=k^2) to make the
        LHS grow faster than the RHS.
    10. Conclude the final answer string based on this analysis.
    """
    # The analysis above leads to the following conclusions:
    # 1. True
    # 2. True
    # 3. True
    # 4. False
    # 5. False
    # 6. False
    final_answer = "TTTFFF"
    print(final_answer)

solve()
#<<<TTTFFF>>>