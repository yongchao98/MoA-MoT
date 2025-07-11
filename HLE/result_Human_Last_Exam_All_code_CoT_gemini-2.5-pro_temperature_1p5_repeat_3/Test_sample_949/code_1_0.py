def solve():
    """
    This function analyzes six statements comparing sums of functions of parenthesis-pair lengths (L) and depths (D).

    My reasoning:
    I established a relationship between L(x), D(x), and properties of a tree representing the parenthesis structure.
    L(x) = 2 * s(x) where s(x) is the number of pairs in the subtree of x.
    D(x) = h(x) where h(x) is the height of the subtree of x.

    To disprove a statement "sum(f(L)) = O(sum(g(D)))", I need a family of strings where the ratio of sums grows infinitely.
    I analyzed extreme cases:
    1. Deeply nested strings (path-like trees): L is proportional to D for all pairs, so the statements hold.
    2. Wide, shallow strings (star-like trees), e.g., s = "(()()...())".
       This creates one pair with large L and small D, and many pairs with small L and D.
       Let's analyze this case:
       - One outer pair: L_out ~ 2n, D_out = 2
       - n inner pairs: L_in = 2, D_in = 1
       For all the given functions f (log, loglog, log^5, L^a with a<1, etc.), f(L_out) grows slower than linear in n.
       The sum on the left (LHS) becomes f(L_out) + n*f(L_in). Because f is sub-linear, n*f(L_in) is the dominant term.
       The sum on the right (RHS) becomes g(D_out) + n*g(D_in), which is also dominated by the term with n.
       Therefore, LHS is O(n) and RHS is O(n). Their ratio is a constant.

    Since these extreme structures do not provide a counterexample, and combinations of them also fail to do so, I conclude that all the given statements are true for all well-formed parenthesis strings.
    """
    answer = "TTTTTT"
    print(answer)

solve()