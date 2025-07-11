import math

def get_sum_of_prime_factors(n, cache={}):
    """
    Calculates the sum of the prime factors of a number n, with multiplicity.
    Uses a cache (dictionary) for memoization to improve performance.
    - For 12 = 2 * 2 * 3, the sum is 2 + 2 + 3 = 7.
    - For 1, the sum is 0.
    """
    if n in cache:
        return cache[n]
    
    if n <= 1:
        return 0

    original_n = n
    total_sum = 0
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            total_sum += d
            temp_n //= d
        d += 1
    if temp_n > 1:
        total_sum += temp_n
    
    cache[original_n] = total_sum
    return total_sum

def find_smallest_sum():
    """
    Finds the smallest integer N that is a sum of two different integers a and b,
    where the sum of the prime factors of a and b is 20.
    """
    target_prime_factor_sum = 20
    
    # Iterate through possible sums N, starting from the smallest (1+2=3).
    # We set a reasonable upper limit, expecting to find the answer quickly.
    for N in range(3, 100): 
        # For each sum N, check all pairs (a, b) where a + b = N and a < b.
        # This inherently satisfies a != b.
        for a in range(1, (N // 2) + 1):
            b = N - a
            
            if a == b:
                continue

            # Check if the prime factor sum condition is met.
            if get_sum_of_prime_factors(a) + get_sum_of_prime_factors(b) == target_prime_factor_sum:
                # We found the smallest N and its corresponding pair (a, b).
                spf_a = get_sum_of_prime_factors(a)
                spf_b = get_sum_of_prime_factors(b)
                
                print(f"Found the smallest number N that satisfies the conditions.")
                print(f"The two integers are a = {a} and b = {b}.")
                print(f"The sum of prime factors for a={a} is {spf_a}.")
                print(f"The sum of prime factors for b={b} is {spf_b}.")
                print(f"The total sum of prime factors is {spf_a} + {spf_b} = {target_prime_factor_sum}.")
                print(f"The final number N is the sum of a and b:")
                print(f"{N} = {a} + {b}")
                print("\nThe required number is:")
                print(N)
                return

    print("No solution found within the checked range.")

if __name__ == '__main__':
    find_smallest_sum()