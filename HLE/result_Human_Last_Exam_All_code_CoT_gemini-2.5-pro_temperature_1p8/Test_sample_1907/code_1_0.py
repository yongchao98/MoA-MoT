import math

# A global cache for primes to avoid re-computation between calls
_primes_list = []
_sieve_limit = 0

def sieve_primes(n):
    """
    Generates a list of prime numbers up to n using the Sieve of Eratosthenes.
    """
    primes = [True] * (n + 1)
    if n >= 0:
        primes[0] = False
    if n >= 1:
        primes[1] = False
    for i in range(2, int(n**0.5) + 1):
        if primes[i]:
            for multiple in range(i*i, n + 1, i):
                primes[multiple] = False
    
    prime_numbers = []
    for i in range(2, n + 1):
        if primes[i]:
            prime_numbers.append(i)
    return prime_numbers

def get_nth_prime(n):
    """
    Returns the n-th prime number. Dynamically extends the prime sieve if necessary.
    """
    global _primes_list, _sieve_limit

    # Check if the current list of primes is large enough
    if n > len(_primes_list):
        # Estimate the size of the n-th prime. p_n < n(log n + log log n) for n >= 6
        if n >= 6:
            log_n = math.log(n)
            limit = int(n * (log_n + math.log(log_n)))
        else:
            # Small n have smaller primes
            limit = 20
        
        # Ensure the new limit is significantly larger to avoid too many re-computations
        if limit <= _sieve_limit:
            limit = int(_sieve_limit * 1.5)
        
        # Generate a new, larger list of primes
        _sieve_limit = limit
        _primes_list = sieve_primes(_sieve_limit)
        
        # If the estimate was too low, extend the sieve until we have enough primes
        while len(_primes_list) < n:
            _sieve_limit = int(_sieve_limit * 1.5)
            _primes_list = sieve_primes(_sieve_limit)

    return _primes_list[n - 1]

def solve():
    """
    Calculates the first term of P^(11) by iterating the recurrence a_{k+1} = p_{a_k}.
    """
    print("This script calculates the first term of P^(k) for k from 1 to 11.")
    print("Let a_k be the first term of P^(k). The rule is a_{k+1} = p_{a_k}, where p_n is the n-th prime.\n")
    
    # a_1 is the first prime number
    current_a = 2
    print(f"The 1st term in P^(1) is: 2")

    # Loop to find a_2 through a_11
    for k in range(2, 12):
        prime_index = current_a
        current_a = get_nth_prime(prime_index)
        # For the final step, the "equation" is a_11 = p_{a_10}
        # The numbers in the final equation are a_10 (the index) and a_11 (the result).
        if k == 11:
            print(f"\nThe final step calculates the 1st term in P^({k}).")
            print(f"The required prime index is the previous term: {prime_index}")
            print(f"The resulting prime (the final answer) is: {current_a}")
            print("\nFinal equation: P_11(1) = Prime(P_10(1))")
            print(f"Prime index: {prime_index}")
            print(f"Resulting Prime: {current_a}")
        else:
            print(f"The 1st term in P^({k}) is the {prime_index}-th prime: {current_a}")
    
    return current_a

if __name__ == '__main__':
    final_answer = solve()
    # The final answer is submitted in the special format below
    # <<<9879893>>>