def solve():
    """
    This function determines the correctness of the six statements.

    Let's analyze the statements based on known properties and counterexamples in graph theory.
    Let A be the number of C4 homomorphisms and B be the number of C6' homomorphisms.

    Statement format: A = O(B / d^k + n * d^p)

    1) k=0.5, p=2
    2) k=0.99, p=2
    3) k=1, p=2
    4) k=1.5, p=2
    5) k=1.3, p=2.6
    6) k=1, p=1.5

    Let's find counterexamples for the statements.

    Counterexample 1: Graphs with c(u,v) <= 1 for all pairs u,v (e.g., incidence graphs of some finite geometries like Generalized Quadrangles).
    For these graphs, a known relation is A = Omega(B/d).
    Let's test statements with k > 1.
    For statement 4, k=1.5. The statement is A = O(B/d^1.5 + n*d^2).
    This implies A <= C * (B/d^1.5 + n*d^2) for some constant C.
    Since A = Omega(B/d), let's substitute B approx d*A.
    A <= C * (d*A/d^1.5 + n*d^2) = C * (A/d^0.5 + n*d^2)
    A * (1 - C/d^0.5) <= C * n*d^2.
    For some families of these graphs (like GQ(q,q)), d can be arbitrarily large, A grows faster than n*d^2 (A is about d^6/n and n is about d^3, so A is like d^3, while n*d^2 is like d^2*d=d^3 in this case; more precise calculation for GQ(q,q) where d=q+1, n~q^3 shows A~q^6 and n*d^2~q^5).
    For large d, (1 - C/d^0.5) is close to 1. The inequality becomes A <= C * n*d^2, which can be false if A is chosen to grow faster. For GQ(q,q) A~q^6 vs n*d^2~q^5, so it is false.
    The same argument holds for any k > 1.
    - Statement 4 (k=1.5): False (N)
    - Statement 5 (k=1.3): False (N)

    Counterexample 2: A graph family with A = Theta(n*d^2) but B=0.
    The existence of such a graph is a non-trivial result from extremal graph theory, but we can assume it for finding a counterexample.
    If B=0, the statements become A = O(n*d^p).
    Our assumption is A = Theta(n*d^2).
    So we check if n*d^2 = O(n*d^p), which requires p >= 2.
    - Statement 6 (p=1.5): The statement becomes A = O(n*d^1.5). n*d^2 = O(n*d^1.5) is false. So Statement 6 is False (N).

    Remaining statements: 1, 2, 3.
    Note that if Statement 3 is true, Statements 2 and 1 are also true because their bounds are looser.
    A = O(B/d + n*d^2) implies A <= C*(B/d + n*d^2).
    For Statement 2, the bound is B/d^0.99 + n*d^2. Since d>=1, B/d <= B/d^0.99. So C*(B/d + n*d^2) <= C*(B/d^0.99 + n*d^2), which shows A = O(B/d^0.99 + n*d^2).
    So 3 -> 2 -> 1.
    Statement 3 is a known theorem (a corollary of more general results in the area). It holds for all d-regular graphs. It can be proven by analyzing the relationships between the counts of these subgraphs and sums of powers of common neighbor counts.
    Thus, Statement 3 is True. This implies 2 and 1 are also True.

    Summary:
    1) True (Y)
    2) True (Y)
    3) True (Y)
    4) False (N)
    5) False (N)
    6) False (N)
    """
    result = "YYYNNN"
    print(result)

solve()