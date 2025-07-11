import math

def solve_and_explain():
    """
    This function explains the solution step-by-step and prints the final answer.
    """

    # Part 1: The condition for a (p,q) pair to be good.
    print("--- Condition for a (p,q) pair to be a 'good' pair ---")
    print("A pair (p,q) is 'good' if we can generate any permutation of {1, ..., n} from any other permutation.")
    print("The allowed operations are swapping two numbers x and y if their difference |x-y| is p or q.")
    
    print("\nThis problem can be modeled using graph theory.")
    print("Let G be a graph with vertices V = {1, 2, ..., n}.")
    print("An edge exists between two vertices i and j if they can be swapped, i.e., |i-j| = p or |i-j| = q.")
    
    print("\nThe ability to generate all permutations (the full symmetric group S_n) is equivalent to this graph G being connected.")
    
    print("\nGiven that 1 < p, q <= n are coprime integers (gcd(p,q) = 1), a known theorem in graph theory states that this graph G is connected if and only if p + q <= n.")
    
    print("\nCONCLUSION (Part 1):")
    print("A pair (p,q) with 1 < p <= n, 1 < q <= n, and gcd(p,q)=1 is a 'good' pair if and only if: p + q <= n")
    print("-" * 50)

    # Part 2: The limit of the probability Pr(n).
    print("\n--- The limit of Pr(n) as n -> infinity ---")
    print("Pr(n) is the probability that a pair (p,q) selected randomly from {2,...,n} x {2,...,n} is a good, coprime pair.")

    print("\n1. Total number of pairs (p,q):")
    print("   p can be chosen from {2, ..., n}, which is (n-1) choices.")
    print("   q can be chosen from {2, ..., n}, which is (n-1) choices.")
    print("   Total pairs for large n is asymptotically n^2.")

    print("\n2. Number of 'good' coprime pairs, N(n):")
    print("   We need to count pairs (p,q) that satisfy:")
    print("   a) 1 < p <= n and 1 < q <= n")
    print("   b) gcd(p,q) = 1 (coprime)")
    print("   c) p + q <= n (the 'good' condition)")

    print("\n3. Asymptotic Analysis:")
    print("   For large n, we can approximate this count by considering the density of points in the p-q plane.")
    print("   - The total pairs (p,q) occupy a square region [2,n] x [2,n], with an area of ~n^2.")
    print("   - The condition p + q <= n confines the valid pairs to a triangular region defined by p>=2, q>=2, p+q<=n. The area of this triangle is ~(1/2)n^2.")
    print("   - The probability that a randomly chosen pair from the square also satisfies p+q <= n is the ratio of their areas: (1/2 * n^2) / n^2 = 1/2.")

    print("\n   - The probability that two integers chosen randomly from a large range are coprime is a famous result: 1/zeta(2) = 6 / pi^2.")
    
    print("\n4. The Limit Calculation:")
    print("   The two events (the geometric constraint and the coprime property) are asymptotically independent. Thus, we can multiply their probabilities to find the limit.")
    print("   lim_{n->inf} Pr(n) = P(p+q <= n) * P(gcd(p,q)=1)")
    print("                      = (1/2) * (6 / pi^2)")
    print("                      = 3 / pi^2")

    print("\nCONCLUSION (Part 2):")
    print("The exact limit of Pr(n) as n approaches infinity is given by the equation:")
    
    numerator = 3
    denominator_name = "pi^2"
    pi_val = math.pi
    denominator_val = pi_val**2
    result = numerator / denominator_val
    
    print(f"  lim Pr(n) = {numerator} / {denominator_name} = {numerator} / ({pi_val:.4f})^2 = {numerator} / {denominator_val:.4f} ≈ {result:.6f}")
    
# Execute the function to print the detailed explanation and the result.
solve_and_explain()