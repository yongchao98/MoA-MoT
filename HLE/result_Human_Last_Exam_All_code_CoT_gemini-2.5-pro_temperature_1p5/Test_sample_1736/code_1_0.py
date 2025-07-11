def sum_prime_factors(n):
    """
    Calculates the sum of the prime factors of a number n, with multiplicity.
    For example, for 12 = 2*2*3, the sum is 2+2+3=7.
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

def find_smallest_number():
    """
    Finds the smallest integer N that is a sum of two different integers a and b,
    such that the sum of their prime factors is 20.
    """
    N = 2
    while True:
        N += 1
        # Iterate through possible values of 'a'. We only need to go up to N/2
        # because b = N - a will cover the other half.
        for a in range(1, N // 2 + 1):
            b = N - a
            
            # The integers a and b must be different.
            # This is naturally handled by iterating a up to N/2,
            # except for the case where N is even and a = N/2, making a=b.
            if a == b:
                continue

            # Check if the sum of prime factors is 20
            spf_a = sum_prime_factors(a)
            spf_b = sum_prime_factors(b)
            if spf_a + spf_b == 20:
                print(f"Found the smallest number N: {N}")
                print(f"The two different integers are a={a} and b={b}.")
                print(f"The final equation satisfying all conditions is:")
                print(f"{a} + {b} = {N}")
                return N

# Run the search and find the number
result = find_smallest_number()