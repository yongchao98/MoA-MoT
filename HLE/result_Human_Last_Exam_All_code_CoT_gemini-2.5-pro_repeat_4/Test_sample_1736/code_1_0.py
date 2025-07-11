import math

def get_sum_of_prime_factors(n):
    """
    Calculates the sum of distinct prime factors of an integer n.
    For example, for n=12 (2^2 * 3), the distinct prime factors are {2, 3},
    and their sum is 2 + 3 = 5.
    """
    if n <= 1:
        return 0
    factors = set()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
    if temp_n > 1:
        factors.add(temp_n)
    return sum(factors)

def solve_number_puzzle():
    """
    Finds the smallest integer N which is the sum of two distinct integers a and b,
    where the sum of the prime factors of a and b is 20.
    """
    # We search for numbers up to a reasonable limit to find the smallest candidates.
    # Numbers like 51 (3*17), 65 (5*13), 91 (7*13) suggest a limit > 100 is safe.
    SEARCH_LIMIT = 150
    TARGET_SPF_SUM = 20
    
    # Store all numbers found for each sum-of-prime-factors (spf) value.
    # This helps find the two smallest numbers if a=b for the minimal choice.
    all_nums_for_spf = {}
    for i in range(1, SEARCH_LIMIT):
        spf = get_sum_of_prime_factors(i)
        if spf not in all_nums_for_spf:
            all_nums_for_spf[spf] = []
        all_nums_for_spf[spf].append(i)

    min_sum = float('inf')
    result_a = -1
    result_b = -1

    # Iterate through possible spf values for 'a'. We only need to go up to half the target sum.
    for spf_a in range(TARGET_SPF_SUM // 2 + 1):
        spf_b = TARGET_SPF_SUM - spf_a
        
        if spf_a in all_nums_for_spf and spf_b in all_nums_for_spf:
            a = all_nums_for_spf[spf_a][0]
            
            # If spf for a and b are the same, we need two different numbers.
            if spf_a == spf_b:
                if len(all_nums_for_spf[spf_a]) > 1:
                    b = all_nums_for_spf[spf_a][1]
                else:
                    # Cannot form a distinct pair, so skip.
                    continue
            else:
                b = all_nums_for_spf[spf_b][0]

            # We have a valid pair (a, b) where a != b. Check if it's the minimal sum.
            current_sum = a + b
            if current_sum < min_sum:
                min_sum = current_sum
                result_a = a
                result_b = b

    # Print the final equation as requested.
    # Sorting ensures the output is consistent, e.g., 3 + 17, not 17 + 3.
    if result_a > result_b:
        result_a, result_b = result_b, result_a
        
    print(f"{result_a} + {result_b} = {min_sum}")

solve_number_puzzle()