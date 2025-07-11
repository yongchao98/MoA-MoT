def solve():
    """
    This function determines the truth value of the six statements.
    Let A be the number of C4 homomorphisms and B be the number of C6' homomorphisms.

    Statement 1: A = O(B/d^0.5 + nd^2)
    Statement 2: A = O(B/d^0.99 + nd^2)
    Statement 3: A = O(B/d + nd^2)
    Statement 4: A = O(B/d^1.5 + nd^2)
    Statement 5: A = O(B/d^1.3 + n d^2.6)
    Statement 6: A = O(B/d + nd^1.5)

    Based on known results in extremal graph theory, for any d-regular graph:
    A = O(B/d + nd^2).
    This means statement (3) is true.

    If statement (3) is true, then since for d >= 1, d^k <= d for k < 1,
    B/d <= B/d^k.
    Therefore, A = O(B/d + nd^2) implies A = O(B/d^k + nd^2) for k < 1.
    This makes statements (1) and (2) true as well.

    To disprove statements (4), (5), and (6), we can use a family of graphs
    that act as counterexamples. Random d-regular graphs (approximated by Ramanujan graphs)
    are a standard choice. For these graphs, for large n and d:
    A is approximately d^4 + nd^2
    B is approximately d^7/n

    Let's check statements 4, 5, 6 with these approximations by finding a relation between n and d
    that shows the O-notation does not hold.

    Statement 6: A = O(B/d + nd^1.5)
    We check: d^4 + nd^2 = O(d^6/n + nd^1.5).
    Let n = d^3. Then LHS is ~d^5. RHS is O(d^3 + d^4.5). d^5 is not O(d^4.5). So statement (6) is False.

    Statement 4: A = O(B/d^1.5 + nd^2)
    We check: d^4 + nd^2 = O(d^5.5/n + nd^2). We only need to check d^4 = O(d^5.5/n).
    This implies n = O(d^1.5). But d-regular graphs exist for any n > d+1.
    If we pick n = d^1.6, for d large enough, d-regular graphs exist.
    With n=d^1.6, the condition n=O(d^1.5) is violated. A more careful check confirms this is a counterexample. So statement (4) is False.

    Statement 5: A = O(B/d^1.3 + nd^2.6)
    We check: d^4 + nd^2 = O(d^5.7/n + nd^2.6).
    The second terms nd^2 vs nd^2.6 are fine.
    We check d^4 = O(d^5.7/n), which implies n=O(d^1.7).
    Let's pick n=d^1.8. The condition is violated. So statement (5) is False.

    Thus, the first three statements are True (Y) and the last three are False (N).
    """
    result = "YYYNNN"
    print(result)

solve()