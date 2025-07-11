def solve():
    """
    This function determines the truth value of six statements about parenthesis strings.
    
    The analysis is as follows:
    1. False. Counterexample: S_n = (()()...()). The sum for L grows as Theta(n), while the sum for D is Theta(1).
    2. False. Counterexample: S_n. The expression log(log(D(x))) is undefined for inner pairs where D(x)=1.
    3. False. Counterexample: S_n. The sum for L grows as Theta(n), while the sum for D is Theta(1).
    4. False. Counterexample: S_n. The condition fails for the outer pair, as L grows with n while D is constant.
    5. True. The exponent on D (0.11) is larger than on L (0.1). This compensates for the L/D gap in all tested string families (wide and deep).
    6. True. The exponent on D (0.5) is larger than on L (0.25), which is sufficient for the relation to hold.
    """
    
    # The final result is a string of 'T' or 'F' for each statement.
    # Based on the analysis, the answers are F, F, F, F, T, T.
    final_answer = "FFFFTT"
    print(final_answer)

solve()