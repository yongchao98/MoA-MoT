import math

def find_special_prime():
    """
    This function finds the largest prime p < 1,000,000 of the form p=4u+1,
    where u=4v+1, and p, u, v are all prime numbers.

    This reflects the logic that would be used to solve the problem on the
    ancient computer, bounded by its maximum integer size of 999999.
    """
    # Define the search limit based on the computer's architecture.
    LIMIT = 1000000

    # Step 1: Generate all prime numbers up to LIMIT using Sieve of Eratosthenes.
    # This is an efficient way to handle primality tests for a large range of numbers.
    is_prime = [True] * LIMIT
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.sqrt(LIMIT)) + 1):
        if is_prime[i]:
            for multiple in range(i * i, LIMIT, i):
                is_prime[multiple] = False

    # Step 2: Determine the search range for v.
    # The relationship is p = 16v + 5.
    # If p < LIMIT, then 16v < LIMIT - 5, so v < (LIMIT - 5) / 16.
    max_v = (LIMIT - 5) // 16

    # Step 3: Iterate downwards from max_v to find the largest triplet.
    # By searching downwards, the first result we find will be the largest.
    for v in range(max_v, 1, -1):
        # Check if v is a prime using our pre-computed list.
        if is_prime[v]:
            u = 4 * v + 1
            # Check if u is within our limit and is prime.
            if u < LIMIT and is_prime[u]:
                p = 4 * u + 1
                # Check if p is within our limit and is prime.
                if p < LIMIT and is_prime[p]:
                    # Found the largest triplet. Print the numbers of the final
                    # equation in the specified format p:u:v and exit.
                    print(f"{p}:{u}:{v}")
                    return

    # This part would be reached if no such prime triplet is found.
    print("No such prime triplet found within the given limit.")


# Execute the function to find and print the result.
find_special_prime()
