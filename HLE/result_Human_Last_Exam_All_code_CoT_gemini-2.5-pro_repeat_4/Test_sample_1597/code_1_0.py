import math

def is_prime(n):
    """Checks if a number is prime."""
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def count_allowed_pairs():
    """
    Calculates the number of allowed pairs (a,b) with 1 <= a, b <= 1000.
    """
    limit = 1000
    
    # Find all primes up to the limit
    primes = [p for p in range(2, limit + 1) if is_prime(p)]
    
    # A dictionary to store the count of powers for each prime
    # key: prime p, value: number of powers p^k <= limit
    prime_power_counts = {}
    
    for p in primes:
        count = 0
        power = p
        while power <= limit:
            count += 1
            # Prevent overflow and unnecessary computation
            if power > limit // p:
                break
            power *= p
        if count > 0:
            prime_power_counts[p] = count
            
    # Calculate the total number of prime powers > 1
    total_prime_powers = sum(prime_power_counts.values())
    
    # Calculate the sum of the squares of these counts
    sum_of_squares = sum(c**2 for c in prime_power_counts.values())
    
    # The number of non-allowed pairs is (sum C_p)^2 - (sum C_p^2)
    num_non_allowed = total_prime_powers**2 - sum_of_squares
    
    # Total pairs
    total_pairs = limit * limit
    
    # The number of allowed pairs
    num_allowed = total_pairs - num_non_allowed
    
    print(f"The total number of pairs is {limit} * {limit} = {total_pairs}.")
    print(f"The number of integers (>1) up to {limit} that are prime powers is {total_prime_powers}.")
    print(f"The sum of squares of counts of powers for each prime is {sum_of_squares}.")
    print(f"The number of non-allowed pairs is {total_prime_powers}^2 - {sum_of_squares} = {total_prime_powers**2} - {sum_of_squares} = {num_non_allowed}.")
    print(f"The number of allowed pairs is {total_pairs} - {num_non_allowed} = {num_allowed}.")

count_allowed_pairs()
