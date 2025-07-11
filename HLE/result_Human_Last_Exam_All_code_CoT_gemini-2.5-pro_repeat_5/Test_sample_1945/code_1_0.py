import math

def solve_good_pairs_problem():
    """
    This function provides the solution to the good pair problem,
    explaining the condition and the probability limit calculation.
    """

    # --- Part 1: Condition for a (p,q) pair to be "good" ---
    condition_explanation = (
        "A pair (p,q) is 'good' if starting from any permutation of {1, 2, ..., n}, "
        "we can reach any other permutation by swapping numbers whose difference is p or q.\n"
        "This is equivalent to the graph with vertices V={1,...,n} and edges between i and j "
        "if |i-j| in {p,q} being connected.\n\n"
        "Given that 1 < p, q <= n and they are coprime, the necessary and "
        "sufficient condition for the graph to be connected is:\n"
        "p + q <= n"
    )
    print("--- Condition for a Good Pair (p,q) ---")
    print(condition_explanation)
    print("-" * 50)

    # --- Part 2: Limit of the probability Pr(n) ---
    limit_explanation = (
        "\nPr(n) is the probability that a pair (p,q) selected randomly from {2,...,n}x{2,...,n} "
        "is coprime and good.\n"
        "1. The total number of pairs (p,q) is (n-1)^2, which is asymptotically n^2.\n"
        "2. The number of 'good, coprime' pairs satisfies p>=2, q>=2, p+q<=n, and gcd(p,q)=1.\n"
        "   These points lie in a triangular region in the p-q plane with an area of ~n^2/2.\n"
        "3. The density of coprime integer pairs is 6/pi^2.\n"
        "4. The number of favorable pairs for large n is approximately (Area) * (Density) = (n^2/2) * (6/pi^2) = 3n^2/pi^2.\n"
        "5. The probability Pr(n) is the ratio of favorable pairs to total pairs.\n"
        "   lim_{n->inf} Pr(n) = lim_{n->inf} (3n^2/pi^2) / n^2"
    )
    print("\n--- Limit of Pr(n) as n -> infinity ---")
    print(limit_explanation)

    numerator = 3
    denominator_symbol = "pi^2"
    denominator_value = math.pi**2
    limit_value = numerator / denominator_value

    print("\nThe final equation for the limit is: lim Pr(n) = " + str(numerator) + " / " + denominator_symbol)
    print("\nAs requested, the numbers in the final equation are:")
    print(f"The numerator is: {numerator}")
    print(f"The denominator is pi-squared, which evaluates to: {denominator_value:.6f}")
    
    print(f"\nThe exact value of the limit is 3/pi^2, which is approximately {limit_value:.6f}.")

solve_good_pairs_problem()