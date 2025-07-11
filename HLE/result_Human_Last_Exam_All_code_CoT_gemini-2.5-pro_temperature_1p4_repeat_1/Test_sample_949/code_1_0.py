def solve():
    """
    Solves the parenthesis problem by analyzing each statement.
    The code prints the reasoning for each statement and then the final answer string.
    """

    # For questions 1, 2, 3, the functions on the RHS (log, loglog, log^5) are 0 or undefined for D=1.
    # We can construct a "comb" string S_k = `(()()...())` with k inner `()` pairs.
    # For this string, there is one root with D=2 and k leaves with D=1.
    # The LHS sum grows with k, while the RHS sum remains constant.
    # Thus, statement_LHS is not O(statement_RHS).
    # So, 1, 2, and 3 are False.
    q1_answer = 'F'
    q2_answer = 'F'
    q3_answer = 'F'

    # For questions 4, 5, 6, the functions on the RHS are positive for D=1.
    # The "comb" string does not work as a counterexample because both sides grow linearly with k.
    # For a "stick" string S_k = `((...()))` with k pairs, for any pair x, L(x) = 2*D(x).
    #
    # Q4: sum(2^sqrt(log(2D))) vs sum(2^sqrt(log(D))). The ratio of terms approaches 1, so the statement holds.
    #
    # Q5: sum((2D)^0.1) vs sum(D^0.11). This is sum(D^0.1) vs sum(D^0.11).
    # The sum grows like integral of x^p, which is x^(p+1).
    # So we compare D^1.1 vs D^1.11. The statement holds.
    #
    # Q6: sum((2D)^0.25) vs sum(D^0.5). This is sum(D^0.25) vs sum(D^0.5).
    # We compare D^1.25 vs D^1.5. The statement holds.
    # It appears that 4, 5, and 6 are True.
    q4_answer = 'T'
    q5_answer = 'T'
    q6_answer = 'T'
    
    final_answer = q1_answer + q2_answer + q3_answer + q4_answer + q5_answer + q6_answer
    print(final_answer)

solve()