def solve():
    """
    This function determines the truth value of six statements about parenthesis strings.
    
    1. sum(log L(x)) = O(sum(log D(x))) -> False. Counterexample: Comb string. LHS grows, RHS is constant.
    2. sum(loglog L(x)) = O(sum(loglog D(x))) -> False. Counterexample: Comb with (()) inside. LHS grows, RHS is constant.
    3. sum(log^5 L(x)) = O(sum(log^5 D(x))) -> False. Counterexample: Comb string. LHS grows, RHS is constant.
    4. sum(2^sqrt(log L(x))) = O(sum(2^sqrt(log D(x)))) -> False. Counterexample: Fat Tree. LHS grows much faster than RHS.
    5. sum(L(x)^0.1) = O(sum(D(x)^0.11)) -> False. Counterexample: Hybrid chain-comb string.
    6. sum(L(x)^(1/4)) = O(sum(D(x)^(1/2))) -> True. The larger exponent for D is sufficient to maintain the inequality for all string structures.
    """
    answer = "FFFFFT"
    print(answer)

solve()