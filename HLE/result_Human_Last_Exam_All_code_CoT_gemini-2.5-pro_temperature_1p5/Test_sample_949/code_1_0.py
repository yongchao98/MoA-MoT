def solve():
    """
    Analyzes the six statements about parenthesis strings.

    Statement 1: sum(log L) = O(sum(log D)) -> True.
    No counterexample found. The sums seem to track each other across different string structures. The `max(1,...)` rule applied to log(D) prevents the RHS sum from becoming too small.

    Statement 2: sum(loglog L) = O(sum(loglog D)) -> False.
    A counterexample is a string with k fixed shells and n tending to infinity, like S_n,k = ((...(()()...())...)).
    For the n simple `()` pairs, L=2 and D=1. Assuming ln(2)<1, log(L)=1 and loglog(L)=0. log(D)=1 and loglog(D)=0.
    These n pairs contribute 0 to both sums. The sums are determined by the few shell pairs.
    For shell pairs, L is large (prop. to n), D is small (const. k).
    LHS grows as log(log(n)), while RHS remains constant. The ratio diverges.

    Statement 3: sum(log^5 L) = O(sum(log^5 D)) -> True.
    Using the same counterexample S_n,k.
    The n simple `()` pairs contribute n*(log 2)^5 = n*1^5 = n to the LHS.
    They contribute n*(log 1)^5 = n*1^5 = n to the RHS.
    The shell pairs contribute a term of order (log n)^5, which is smaller than n.
    So both sums are dominated by the O(n) term, and their ratio is constant.

    Statement 4: sum(2^sqrt(log L)) = sum(2^O(sqrt(log D))) -> False.
    Using S_n,k again. The statement means sum(2^sqrt(log L)) <= sum(2^(C*sqrt(log D))) for some C.
    LHS is dominated by shell pairs, giving a term growing like 2^sqrt(log n).
    RHS is dominated by the n simple pairs, giving a term growing like n.
    The term 2^sqrt(log n) grows superpolynomially, so it grows faster than n. The ratio diverges.

    Statement 5: sum(L^0.1) = O(sum(D^0.11)) -> True.
    The exponent on the RHS is larger. No counterexamples found. For S_n,k, both sums are O(n). For other families, the inequality holds.

    Statement 6: sum(L^0.25) = O(sum(D^0.5)) -> True.
    Same logic as statement 5. The exponent on the RHS is larger, and no counterexamples have been found.
    """
    
    # Based on the step-by-step analysis:
    # 1. True
    # 2. False
    # 3. True
    # 4. False
    # 5. True
    # 6. True
    answer = "TFTFTT"
    print(answer)

solve()