import math

class PrimeChecker:
    """
    A data structure for efficient prime number operations using a pre-computed sieve.
    - isprime(p): O(1)
    - primes(n): O(n)
    """
    def __init__(self, max_num=10000):
        self._max_num = max_num
        # In a real memory-optimized implementation (e.g., in C or using a library),
        # this would be a bit array. For Python, a list of booleans is more readable.
        self._sieve = [True] * self._max_num
        self._sieve[0] = self._sieve[1] = False

        # Sieve of Eratosthenes to pre-compute all primes
        for i in range(2, int(math.sqrt(self._max_num))):
            if self._sieve[i]:
                # Mark all multiples of i as not prime
                for multiple in range(i*i, self._max_num, i):
                    self._sieve[multiple] = False

    def isprime(self, p):
        """
        Checks if p is a prime number in O(1) time.
        p must be less than the initial max_num.
        """
        if 0 <= p < self._max_num:
            return self._sieve[p]
        raise ValueError(f"Input {p} is out of the pre-computed range [0, {self._max_num-1}]")

    def primes(self, n):
        """
        Returns a list of all primes p <= n in O(n) time.
        n must be less than the initial max_num.
        """
        if 0 <= n < self._max_num:
            return [i for i in range(n + 1) if self._sieve[i]]
        raise ValueError(f"Input {n} is out of the pre-computed range [0, {self._max_num-1}]")

    def get_data_structure_size_in_bytes(self):
        """
        Calculates the size of the most memory-efficient implementation (a bit array).
        """
        # We need one bit for each number from 0 to max_num - 1.
        num_bits = self._max_num
        # There are 8 bits in a byte. Use ceiling division.
        num_bytes = math.ceil(num_bits / 8)
        return int(num_bytes)

# Main execution block
if __name__ == "__main__":
    # The maximum value for p and n is less than 10000.
    MAX_LIMIT = 10000
    
    # Initialize the data structure
    prime_checker = PrimeChecker(MAX_LIMIT)

    # Example Usage:
    # 1. isprime(p)
    p = 9973
    print(f"Is {p} a prime? {prime_checker.isprime(p)}")
    
    p = 9975
    print(f"Is {p} a prime? {prime_checker.isprime(p)}")

    # 2. primes(n)
    n = 30
    print(f"Primes up to {n}: {prime_checker.primes(n)}")

    # Calculate and print the maximal size of the data structure
    size = prime_checker.get_data_structure_size_in_bytes()
    
    print("\n--- Data Structure Size Calculation ---")
    print(f"To store primality for {MAX_LIMIT} numbers, we need {MAX_LIMIT} bits.")
    print("Since there are 8 bits in a byte, the required memory is:")
    print(f"ceiling({MAX_LIMIT} / 8) = {size} bytes.")
    print(f"The maximal size of this data structure is {size} bytes.")

<<<1250>>>