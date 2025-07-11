import math

memo_spf = {}

def get_sum_of_prime_factors(n):
    """
    Calculates the sum of the unique prime factors of an integer n.
    Results are cached in a dictionary to improve performance.
    """
    if n in memo_spf:
        return memo_spf[n]
    
    if n <= 1:
        return 0

    original_n = n
    factors = set()
    d = 2
    
    while d * d <= n:
        if n % d == 0:
            factors.add(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.add(n)

    result = sum(factors)
    memo_spf[original_n] = result
    return result

def find_the_number():
    """
    Searches for the smallest integer N = a + b, where a != b,
    and the sum of the prime factors of a and b is 20.
    """
    min_n = float('inf')
    result_a = -1
    result_b = -1

    # A search limit of 100 for a and b is sufficient, as the smallest N
    # will likely be formed from relatively small numbers.
    search_limit = 100

    for a in range(2, search_limit):
        for b in range(a + 1, search_limit):
            spf_a = get_sum_of_prime_factors(a)
            spf_b = get_sum_of_prime_factors(b)

            if spf_a + spf_b == 20:
                current_n = a + b
                if current_n < min_n:
                    min_n = current_n
                    result_a = a
                    result_b = b

    if result_a != -1:
        print(f"The number that satisfies the conditions is N = {min_n}.")
        print("This number is the result of the following equation:")
        # Output each number in the final equation as requested.
        print(f"{result_a} + {result_b} = {min_n}")
    else:
        print("A solution could not be found within the search limit.")

if __name__ == '__main__':
    find_the_number()