import math

class PrimeSieve:
    """
    A data structure to efficiently check for primality and list primes
    up to a pre-computed limit.
    """
    def __init__(self, max_limit):
        """
        Initializes the sieve by pre-computing primes up to max_limit.
        """
        self._max_limit = max_limit
        # Create a boolean list to store primality. True means potentially prime.
        self._sieve = [True] * (max_limit + 1)
        if max_limit >= 0:
            self._sieve[0] = False
        if max_limit >= 1:
            self._sieve[1] = False
        
        # Sieve of Eratosthenes algorithm
        for i in range(2, int(math.sqrt(max_limit)) + 1):
            if self._sieve[i]:
                # Mark all multiples of i as not prime
                for multiple in range(i * i, max_limit + 1, i):
                    self._sieve[multiple] = False

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be within the pre-computed limit.
        """
        if not 0 <= p <= self._max_limit:
            raise ValueError(f"Input {p} is out of the pre-computed range [0, {self._max_limit}]")
        return self._sieve[p]

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n in O(n) time.
        n must be within the pre-computed limit.
        """
        if not 0 <= n <= self._max_limit:
            raise ValueError(f"Input {n} is out of the pre-computed range [0, {self._max_limit}]")
        
        prime_list = []
        for i in range(n + 1):
            if self._sieve[i]:
                prime_list.append(i)
        return prime_list

# Main execution
if __name__ == "__main__":
    # The problem specifies a limit of less than 10000, so we compute up to 9999.
    MAX_N = 9999
    
    # Instantiate the data structure
    prime_checker = PrimeSieve(MAX_N)

    # --- Demonstrate the operations ---
    # 1. isprime(p)
    p = 9973 # A large prime
    print(f"Is {p} a prime number? {prime_checker.isprime(p)}")
    p = 9975
    print(f"Is {p} a prime number? {prime_checker.isprime(p)}")
    
    # 2. primes(n)
    n = 30
    print(f"\nPrimes less than or equal to {n}: {prime_checker.primes(n)}")

    # --- Calculate and print the maximal size ---
    print("\n--- Maximal Size Calculation ---")
    print("The most memory-efficient design uses a bit array, where each number")
    print("from 0 to 9999 is represented by a single bit.")
    
    total_numbers = 10000  # For numbers 0 to 9999
    bits_per_byte = 8
    size_in_bytes = total_numbers // bits_per_byte
    
    print("\nThe final equation for the size in bytes is:")
    print(f"{total_numbers} / {bits_per_byte} = {size_in_bytes}")

<<<1250>>>