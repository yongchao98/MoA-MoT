def solve():
    """
    Solves the parenthesis string property problem.

    This function analyzes six statements about the length (L) and depth (D)
    of matched parenthesis pairs in a string. It determines whether each statement
    is True for all strings (T) or False (F).

    The analysis leads to the following conclusions:
    1. False. A 'wide' string like `(()()...())` serves as a counterexample.
    2. False. The same 'wide' string is a counterexample, considering the domains of the log-log function.
    3. False. The same 'wide' string is a counterexample.
    4. True. The growth of the function is slow enough that the sums are dominated by the number of pairs, leading to comparable growth for both sides.
    5. True. The exponent on the D-side is larger than on the L-side, which is sufficient to maintain the big-O relationship for all string structures.
    6. True. Similar to statement 5, the exponent on the D-side (0.5) is larger than on the L-side (0.25).
    """
    
    # The final result is a string concatenation of T/F findings.
    answer = "FFFTTT"
    print(answer)

solve()