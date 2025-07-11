import sys

class PrimeDataStructure:
    """
    A data structure to efficiently check for primality and list primes
    up to a certain limit.
    """
    def __init__(self, max_num=9999):
        """
        Initializes the data structure by pre-computing primes using a sieve.
        The data structure used is a boolean list, representing a sieve.
        A bit array would be more memory-efficient but a list is used for
        readability.
        """
        self._max_num = max_num
        # Create a boolean list to store primality information.
        # self._sieve[i] will be True if i is prime, False otherwise.
        self._sieve = [True] * (max_num + 1)
        if max_num >= 0:
            self._sieve[0] = False
        if max_num >= 1:
            self._sieve[1] = False

        # Sieve of Eratosthenes algorithm
        for p in range(2, int(max_num**0.5) + 1):
            if self._sieve[p]:
                # Mark all multiples of p as not prime
                for i in range(p * p, max_num + 1, p):
                    self._sieve[i] = False

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        """
        if 0 <= p <= self._max_num:
            return self._sieve[p]
        return False

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n in O(n) time.
        """
        if n > self._max_num:
            n = self._max_num
        
        prime_list = []
        for i in range(2, n + 1):
            if self._sieve[i]:
                prime_list.append(i)
        return prime_list

if __name__ == "__main__":
    # The maximum number for which we need to store primality is 9999.
    # Therefore, we need to store 10000 boolean values (for numbers 0 through 9999).
    prime_checker = PrimeDataStructure(max_num=9999)

    # Example Usage:
    print(f"Is 9973 a prime number? {prime_checker.isprime(9973)}")
    print(f"Is 100 a prime number? {prime_checker.isprime(100)}")
    print(f"Primes up to 30: {prime_checker.primes(30)}")
    print("-" * 20)

    # The most memory-efficient design uses a bit array.
    # We calculate the size of this efficient data structure in bytes.
    print("Calculating the size of the most memory-efficient data structure (bit array):")
    
    # We need to store information for numbers 0 through 9999.
    numbers_to_store = 10000
    bits_per_byte = 8
    
    # Each number's primality can be stored in a single bit.
    size_in_bytes = numbers_to_store / bits_per_byte
    
    # The problem asks to output each number in the final equation.
    print("Final equation for size calculation:")
    print(f"{numbers_to_store} (bits) / {bits_per_byte} (bits per byte) = {int(size_in_bytes)} (bytes)")

<<<1250>>>