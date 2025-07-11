def solve_and_explain():
    """
    This function prints the derivation and the final answer to the problem.
    """
    print("--- Part 1: Condition for a Good Pair (p,q) ---")
    print("A pair (p,q) is 'good' if the graph on vertices {1,...,n} with edges connecting numbers with difference p or q is connected.")
    print("We are given that p and q are coprime, so gcd(p,q) = 1.")
    print("A known theorem states that this graph is connected if and only if n >= p + q - gcd(p,q).")
    print("Substituting gcd(p,q) = 1, the condition becomes:")
    print("n >= p + q - 1, which is equivalent to p + q <= n + 1.")
    print("\nCondition: A coprime pair (p,q) with 1 < p,q <= n is good if and only if p + q <= n + 1.\n")

    print("--- Part 2: The Limit of the Probability Pr(n) ---")
    print("Pr(n) is the probability that a randomly chosen coprime pair (p,q) with 1 < p,q <= n is a good pair.")
    print("Let C_n be the set of coprime pairs (p,q) with 1 < p,q <= n.")
    print("Let G_n be the set of good pairs, which are pairs in C_n satisfying p + q <= n + 1.")
    print("The probability is Pr(n) = |G_n| / |C_n|.\n")
    print("To find the limit as n -> infinity, we compare the areas of the regions where these points lie.")
    print("The density of coprime pairs in a large region of the integer grid is uniform (6/pi^2).")
    print("So, lim Pr(n) = lim Area(Region for G_n) / Area(Region for C_n).\n")
    print("Region for C_n: A square defined by 2 <= p <= n and 2 <= q <= n.")
    print("Area(C_n) = (n-1) * (n-1) = (n-1)^2.\n")
    print("Region for G_n: A triangle defined by 2 <= p, 2 <= q, and p + q <= n + 1.")
    print("Area(G_n) = 1/2 * base * height = 1/2 * (n-3) * (n-3) = 1/2 * (n-3)^2.\n")
    
    print("The limit calculation is as follows:")
    print("lim_{n->inf} Pr(n) = lim_{n->inf} Area(G_n) / Area(C_n)")
    print("                 = lim_{n->inf} [1/2 * (n-3)^2] / [(n-1)^2]")
    print("                 = 1/2 * lim_{n->inf} (n^2 - 6*n + 9) / (n^2 - 2*n + 1)")
    print("                 = 1/2 * lim_{n->inf} (1 - 6/n + 9/n^2) / (1 - 2/n + 1/n^2)")
    print("                 = 1/2 * (1 / 1)")
    print("                 = 0.5\n")

    limit_value = 0.5
    print(f"The exact value of lim_{{n->inf}} Pr(n) is {limit_value}.")

if __name__ == '__main__':
    solve_and_explain()