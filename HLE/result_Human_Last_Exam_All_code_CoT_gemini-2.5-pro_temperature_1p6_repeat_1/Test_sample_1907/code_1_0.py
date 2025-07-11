import math

class PrimeHandler:
    """
    A class to handle prime number calculations efficiently.
    It uses a dynamic sieve that extends as needed.
    """
    def __init__(self):
        # Initialize with the first few primes to bootstrap the process
        self._primes = [2, 3, 5, 7, 11, 13]
        self._sieve_limit = 13

    def _sieve_upto(self, limit):
        """
        Generates primes up to a given limit using a sieve and stores them.
        """
        if limit <= self._sieve_limit:
            return
        
        self._sieve_limit = limit
        is_prime = [True] * (limit + 1)
        is_prime[0] = is_prime[1] = False
        
        for i in range(2, int(math.sqrt(limit)) + 1):
            if is_prime[i]:
                for multiple in range(i*i, limit + 1, i):
                    is_prime[multiple] = False
        
        self._primes = [p for p, is_p in enumerate(is_prime) if is_p]

    def get_nth_prime(self, n):
        """
        Gets the n-th prime number (1-indexed).
        Estimates the required sieve size and extends it if necessary.
        """
        if n <= 0:
            raise ValueError("n must be a positive integer")
        
        while n > len(self._primes):
            # Estimate the upper bound for the n-th prime.
            # Rosser's theorem states p_n < n * (ln(n) + ln(ln(n))) for n >= 6.
            # We add a small margin to be safe.
            if n < 6:
                # Use a generous default for small n where approximation is poor.
                limit = 30
            else:
                log_n = math.log(n)
                log_log_n = math.log(log_n)
                limit_approx = n * (log_n + log_log_n)
                limit = int(limit_approx * 1.15) 

            # Extend sieve to at least double the current size to avoid frequent re-sieving.
            new_limit = max(self._sieve_limit * 2, limit)
            self._sieve_upto(new_limit)

        return self._primes[n-1]

def solve():
    """
    Solves the problem by calculating the sequence a_k and printing the results.
    """
    print("Let p(n) be the n-th prime number.")
    print("Let a(k) be the 1st term in P^(k).")
    print("The problem definition leads to the recurrence relation: a(k+1) = p(a(k)).")
    print("The 1st term of P^(1) (the set of primes) is a(1) = 2.")
    print("We need to find a(11).\n")

    prime_handler = PrimeHandler()
    
    # Store the sequence a_1, a_2, ..., a_11
    values = []
    a = 2
    values.append(a)

    for _ in range(10):
        a = prime_handler.get_nth_prime(a)
        values.append(a)

    print("The calculation steps are:")
    print(f"a(1) = 2")
    for k in range(2, 12):
        print(f"a({k}) = p(a({k-1})) = p({values[k-2]}) = {values[k-1]}")

    print("\nThe full equation showing each number substitution is:")
    
    parts = []
    num_p_nestings = 10
    for i in range(10): # Corresponds to a_1 through a_10
        prefix = "p(" * num_p_nestings
        suffix = ")" * num_p_nestings
        part = f"{prefix}{values[i]}{suffix}"
        parts.append(part)
        num_p_nestings -= 1
        
    final_result = values[10]
    full_equation = " = ".join(parts)
    full_equation += f" = {final_result}"

    print(full_equation)
    print(f"\nThus, the 1st term in P^(11) is {final_result}.")


if __name__ == "__main__":
    solve()
    print("\n<<<10166441>>>")
