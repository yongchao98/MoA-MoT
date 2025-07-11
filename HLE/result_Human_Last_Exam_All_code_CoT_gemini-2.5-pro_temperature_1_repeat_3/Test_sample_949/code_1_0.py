def solve():
    """
    This function determines the truth value of six statements about parenthesis strings.
    My reasoning is as follows:
    I test each statement against a specific family of strings: S_k = "(()()...())" with k inner pairs.
    This string has one outer pair (x_0) with L=2(k+1), D=2, and k inner pairs (x_i) with L=2, D=1.

    1. sum(log L) vs sum(log D): LHS is O(k), RHS is O(1) because log(1)=0. False.
    2. sum(loglog L) vs sum(loglog D): LHS is O(k), RHS is O(1) because loglog(1) is undefined, so the term is omitted. False.
    3. sum(log^5 L) vs sum(log^5 D): LHS is O(k), RHS is O(1) because log^5(1)=0. False.
    4. sum(2^sqrt(log L)) vs sum(2^sqrt(log D)): LHS is O(k), RHS is O(k) because 2^sqrt(log 1) = 1. The counterexample fails. True.
    5. sum(L^0.1) vs sum(D^0.11): LHS is O(k), RHS is O(k) because 1^0.11 = 1. The counterexample fails. True.
    6. sum(L^0.25) vs sum(D^0.5): LHS is O(k), RHS is O(k) because 1^0.5 = 1. The counterexample fails. True.
    
    The final answer is a string concatenating the truth values (T/F).
    """
    # The reasoning provided above leads to the following conclusions for the statements:
    # 1. False
    # 2. False
    # 3. False
    # 4. True
    # 5. True
    # 6. True
    answer = "FFFTTT"
    print(answer)

solve()