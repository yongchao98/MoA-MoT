def solve():
    """
    This function determines the truth value of six statements about parenthesis strings.
    
    My reasoning is as follows:
    1. False. Consider S_k = "(()()...())" with k inner pairs. Sum(log L) ~ k, while Sum(log D) is constant.
    2. False. For S_k, Sum(log log L) ~ k, while Sum(log log D) is constant (sum over D>1).
    3. False. For S_k, Sum(log^5 L) ~ (log k)^5, while Sum(log^5 D) is constant.
    4. False. For S_k, Sum(2^sqrt(log L)) ~ k, while Sum(2^sqrt(log D)) is constant.
    5. True. The exponent on D (0.11) is larger than on L (0.1). This advantage is sufficient for the sum over D to keep pace with the sum over L, even for extreme string structures. For wide strings, both sums are O(k). For deep strings, Sum(L^0.1) ~ k^1.1 and Sum(D^0.11) ~ k^1.11.
    6. True. The exponent on D (0.5) is much larger than on L (0.25), making this relationship hold even more strongly than in statement 5.
    """
    answer = "FFFFTT"
    print(answer)

solve()