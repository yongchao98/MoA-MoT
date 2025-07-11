def solve():
    """
    Analyzes six statements comparing sums over lengths and depths of parenthesis pairs.
    
    Let's denote the statements as S1 to S6.
    
    A key insight is the relationship L(x) = 2 * (number of pairs in the subtree rooted at x)
    and D(x) = (height of the subtree rooted at x).
    Let n_x be the number of pairs inside x. Then L(x) = 2 * (1 + n_x).
    Also, D(x) <= 1 + n_x.
    This implies L(x) >= 2D(x).

    We analyze by considering two main families of strings (trees):
    1. Path graphs (deeply nested): s = ((...())), L(x) is proportional to D(x).
    2. Star graphs (wide and shallow): s = (()()...()), L(x) can be much larger than D(x) for the root, but for most pairs (leaves), L(x) and D(x) are small constants.

    S1: sum(log L) vs sum(log D)
    - Path: sum(log(2i)) vs sum(log i). Both are Theta(N log N). Holds.
    - Star: sum is dominated by N leaf nodes. LHS is ~ N*log(2), RHS is ~ N*log(1) (clarified to N*1). Both are Theta(N). Holds.
    - Conclusion: True.

    S2: sum(loglog L) vs sum(loglog D)
    - The loglog function grows extremely slowly. The logic is similar to S1.
    - Path: sum(loglog(2i)) vs sum(loglog i). Both are Theta(N loglog N). Holds.
    - Star: The sum is dominated by the N leaf terms. With clarification, both LHS and RHS are Theta(N). Holds.
    - Conclusion: True.

    S3: sum(log^5 L) vs sum(log^5 D)
    - Similar to S1.
    - Path: sum((log 2i)^5) vs sum((log i)^5). Both are Theta(N (log N)^5). Holds.
    - Star: Both sums are Theta(N). Holds.
    - Conclusion: True.

    S4: sum(2^sqrt(log L)) vs sum(2^O(sqrt(log D)))
    - The function 2^sqrt(log t) grows slower than any power t^epsilon.
    - Path: sum(2^sqrt(log 2i)) vs sum(2^sqrt(log i)). Asymptotically similar. Holds.
    - Star: Both sums are Theta(N). Holds.
    - Conclusion: True.

    S5: sum(L^0.1) vs sum(D^0.11)
    - Here we compare powers p=0.1 and q=0.11.
    - Path: sum((2i)^0.1) vs sum(i^0.11). This is Theta(N^1.1) vs Theta(N^1.11). Since 1.1 < 1.11, it holds.
    - Star: Both sums are Theta(N). Holds.
    - Conclusion: True, crucially because p < q.

    S6: sum(L^0.25) vs sum(D^0.5)
    - Here we compare powers p=0.25 and q=0.5.
    - Path: sum((2i)^0.25) vs sum(i^0.5). This is Theta(N^1.25) vs Theta(N^1.5). Since 1.25 < 1.5, it holds.
    - Star: Both sums are Theta(N). Holds.
    - Conclusion: True, again because p < q.
    
    All statements are true.
    """
    answer = "TTTTTT"
    print(answer)

solve()