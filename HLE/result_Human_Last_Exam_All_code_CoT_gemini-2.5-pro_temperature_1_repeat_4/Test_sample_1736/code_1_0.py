def get_sum_of_prime_factors(n):
    """
    Calculates the sum of the prime factors of a number, including multiplicity.
    For example, for 12 = 2 * 2 * 3, the sum is 2 + 2 + 3 = 7.
    """
    if n <= 1:
        return 0
    
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
    return total_sum

def find_smallest_n():
    """
    Finds the smallest integer N that is the sum of two different integers a and b,
    such that the sum of the prime factors of a and b is 20.
    """
    # Iterate through possible sums N, starting from 3 (1+2).
    # The upper limit of 100 is chosen as a safe bound; the answer is expected to be smaller.
    for n_sum in range(3, 101):
        # Iterate through possible values for 'a'.
        # We only need to check a up to n_sum / 2 to find unique pairs (a, b) where a < b.
        for a in range(1, n_sum // 2 + 1):
            b = n_sum - a
            
            # Condition ¬R: a and b must be different.
            # The loop structure with a < b ensures this, but we explicitly check a != b
            # to handle the case where n_sum is even and a = n_sum / 2.
            if a == b:
                continue

            # Condition ¬Q: Sum of prime factors of a and b is 20.
            sum_factors_a = get_sum_of_prime_factors(a)
            sum_factors_b = get_sum_of_prime_factors(b)
            
            if sum_factors_a + sum_factors_b == 20:
                # Condition ¬P: We found the smallest N because we are iterating n_sum upwards.
                # Print the equation with each number.
                print(f"{a} + {b} = {n_sum}")
                return n_sum

# Execute the search and print the result.
final_number = find_smallest_n()
# The final answer is the number N itself.
print(f"<<<{final_number}>>>")
