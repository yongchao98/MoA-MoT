import math

def count_allowed_pairs(limit):
    """
    Calculates the number of allowed ordered pairs (a,b) with 1 <= a,b <= limit.

    An ordered pair (a,b) is allowed if it is NOT of the form (p^k, q^l) where
    p and q are distinct primes and k, l >= 1.
    """
    
    # Step 1: Generate primes up to the limit using a sieve
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(math.sqrt(limit)) + 1):
        if is_prime[p]:
            for multiple in range(p * p, limit + 1, p):
                is_prime[multiple] = False
    
    primes = [p for p, is_p in enumerate(is_prime) if is_p]

    # Step 2: For each prime, count its powers <= limit
    # c_p in the explanation corresponds to counts_per_prime[p]
    counts_per_prime = {}
    for p in primes:
        count = 0
        power_val = p
        while power_val <= limit:
            count += 1
            # Check for potential overflow before multiplication, although not strictly necessary for limit=1000
            if limit // p < power_val:
                break
            power_val *= p
        if count > 0:
            counts_per_prime[p] = count

    # Step 3: Calculate S (sum of counts) and S2 (sum of squares of counts)
    counts = list(counts_per_prime.values())
    sum_of_counts = sum(counts)
    sum_of_squares = sum(c * c for c in counts)

    # Step 4: Calculate the number of disallowed pairs
    # Disallowed pairs = (sum of counts)^2 - (sum of squares of counts)
    num_disallowed_pairs = sum_of_counts * sum_of_counts - sum_of_squares

    # Step 5: Calculate total pairs and subtract disallowed pairs
    total_pairs = limit * limit
    num_allowed_pairs = total_pairs - num_disallowed_pairs

    print(f"For 1 <= a,b <= {limit}:")
    print(f"Total number of pairs = {limit} * {limit} = {total_pairs}")
    print(f"Number of integers > 1 that are prime powers = {sum_of_counts}")
    print(f"Sum of squares of these counts (sum_p c_p^2) = {sum_of_squares}")
    print(f"Number of disallowed pairs = {sum_of_counts}^2 - {sum_of_squares} = {num_disallowed_pairs}")
    print(f"Number of allowed pairs = {total_pairs} - {num_disallowed_pairs} = {num_allowed_pairs}")

limit = 1000
count_allowed_pairs(limit)

# Final answer extraction
# The question only needs the final numerical answer as per the problem description.
total_pairs = limit * limit
is_prime = [True] * (limit + 1)
is_prime[0] = is_prime[1] = False
for p in range(2, int(math.sqrt(limit)) + 1):
    if is_prime[p]:
        for multiple in range(p * p, limit + 1, p):
            is_prime[multiple] = False
primes = [p for p, is_p in enumerate(is_prime) if is_p]
counts_per_prime = {}
for p in primes:
    count = 0
    power_val = p
    while power_val <= limit:
        count += 1
        if limit // p < power_val:
            break
        power_val *= p
    if count > 0:
        counts_per_prime[p] = count
counts = list(counts_per_prime.values())
sum_of_counts = sum(counts)
sum_of_squares = sum(c * c for c in counts)
num_disallowed_pairs = sum_of_counts * sum_of_counts - sum_of_squares
num_allowed_pairs = total_pairs - num_disallowed_pairs
print(f"\nFinal answer: {num_allowed_pairs}")