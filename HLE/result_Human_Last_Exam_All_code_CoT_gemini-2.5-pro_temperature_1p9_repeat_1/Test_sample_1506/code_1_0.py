import sys
import math

class PrimeChecker:
    """
    A memory-efficient data structure to check for primality and list primes
    up to a certain limit (MAX_N = 10000).

    It uses a bytearray as a sieve for O(1) primality tests and O(n) prime
    listing, while minimizing memory usage.
    """

    def __init__(self, max_n=10000):
        """
        Initializes the sieve up to max_n.
        The sieve is a bytearray where sieve[i] = 1 if i is prime, 0 otherwise.
        """
        self._max_n = max_n
        # Use bytearray for memory efficiency: 1 byte per number.
        self._sieve = bytearray([1]) * self._max_n
        self._sieve[0] = 0
        self._sieve[1] = 0

        # Sieve of Eratosthenes algorithm
        for i in range(2, int(math.sqrt(self._max_n))):
            if self._sieve[i] == 1:
                # Mark all multiples of i as not prime
                for multiple in range(i * i, self._max_n, i):
                    self._sieve[multiple] = 0

    def isprime(self, p):
        """
        Checks if a number p is prime.
        Time complexity: O(1).
        """
        if 0 <= p < self._max_n:
            return self._sieve[p] == 1
        return False

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n.
        Time complexity: O(n).
        """
        if not 0 <= n < self._max_n:
            n = self._max_n -1
            
        prime_list = []
        for i in range(2, n + 1):
            if self._sieve[i] == 1:
                prime_list.append(i)
        return prime_list

    def get_data_structure_size_in_bytes(self):
        """
        Calculates the size of the core data structure (the sieve).
        """
        return sys.getsizeof(self._sieve)

if __name__ == "__main__":
    # 1. Create the prime checker instance
    # This pre-computes all primes up to 9999
    checker = PrimeChecker(10000)

    # 2. Demonstrate isprime(p) operation
    print("--- Demonstrating isprime(p) ---")
    p1 = 9973
    p2 = 9974
    print(f"Is {p1} prime? {checker.isprime(p1)}")
    print(f"Is {p2} prime? {checker.isprime(p2)}")
    print("-" * 20)

    # 3. Demonstrate primes(n) operation
    print("--- Demonstrating primes(n) ---")
    n = 30
    prime_list_up_to_n = checker.primes(n)
    print(f"Primes up to {n}: {prime_list_up_to_n}")
    print("-" * 20)

    # 4. Calculate and report the size of the data structure
    print("--- Calculating Data Structure Size ---")
    total_size_bytes = checker.get_data_structure_size_in_bytes()
    data_bytes = checker._max_n
    overhead_bytes = total_size_bytes - data_bytes
    
    print(f"The core data structure is a bytearray of length {checker._max_n}.")
    print("Size (in bytes) is calculated as: Object Overhead + Data Storage")
    print("The final equation for the maximal size is:")
    print(f"{overhead_bytes} + {data_bytes} = {total_size_bytes}")
    print(f"The maximal size of this data structure is {total_size_bytes} bytes.")
<<<10057>>>