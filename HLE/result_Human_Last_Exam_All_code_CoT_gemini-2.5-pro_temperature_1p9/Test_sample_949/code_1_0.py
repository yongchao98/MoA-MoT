def solve():
    """
    Solves the problem by determining the truth value of each statement.

    1. False. Consider a tree of depth d=1 with k leaves. As k->inf, LHS/RHS ratio is bounded.
       Consider a tree of arity k and depth d. As k->inf, ratio is O(log k). False.
    2. False. As k->inf, ratio is O(log log k). False.
    3. False. As k->inf, ratio is O((log k)^5). False.
    4. False. Consider a tree of arity k and depth d. L_root ~ k^d, D_root = d.
       Term for root on LHS is 2^sqrt(d log k), on RHS is 2^sqrt(log d). The ratio blows up. False.
    5. True. p=0.1, q=0.11. q>p. Analysis on various string structures (deep, wide, k-ary trees)
       shows the ratio of sums is always bounded.
    6. True. p=0.25, q=0.5. q>p. Same reasoning as 5.
    """
    answer = "FFFFTT"
    print(answer)

solve()