import math

class PrimeDataStructure:
    """
    A data structure to efficiently check for primality and list primes
    up to a pre-computed limit.
    """

    def __init__(self, max_limit=10000):
        """
        Initializes the data structure by pre-computing primes using the Sieve of Eratosthenes.
        """
        self._max_limit = max_limit
        # In a real memory-optimized implementation, this would be a bit array.
        # For clarity in Python, we use a list of booleans.
        self._sieve = [True] * self._max_limit
        self._sieve[0] = self._sieve[1] = False
        for i in range(2, int(math.sqrt(self._max_limit))):
            if self._sieve[i]:
                # Mark all multiples of i as not prime
                for multiple in range(i*i, self._max_limit, i):
                    self._sieve[multiple] = False

    def isprime(self, p):
        """
        Checks if p is a prime number in O(1) time.
        p must be less than the max_limit.
        """
        if not 0 <= p < self._max_limit:
            raise ValueError(f"Input p must be between 0 and {self._max_limit-1}")
        return self._sieve[p]

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n in O(n) time.
        n must be less than the max_limit.
        """
        if not 0 <= n < self._max_limit:
            raise ValueError(f"Input n must be between 0 and {self._max_limit-1}")
        
        prime_list = []
        for i in range(n + 1):
            if self._sieve[i]:
                prime_list.append(i)
        return prime_list

def calculate_max_size():
    """
    Calculates the maximal size of the most efficient data structure (a bit array).
    """
    numbers_to_store = 10000  # For numbers 0 to 9999
    bits_per_number = 1
    total_bits = numbers_to_store * bits_per_number
    # There are 8 bits in a byte
    total_bytes = total_bits / 8
    
    print("Calculation for the maximal size of the data structure:")
    print(f"Numbers to store: {numbers_to_store} (for 0 to 9999)")
    print(f"Bits per number: {bits_per_number}")
    print(f"Total bits required: {total_bits}")
    print(f"Total bytes required: {total_bits} / 8 = {int(total_bytes)}")
    print(f"The maximal size of this data structure is {int(total_bytes)} bytes.")


if __name__ == '__main__':
    # Demonstrate the data structure
    prime_checker = PrimeDataStructure()

    # --- isprime(p) demonstration ---
    p1 = 29
    p2 = 30
    print(f"Is {p1} prime? {prime_checker.isprime(p1)}")
    print(f"Is {p2} prime? {prime_checker.isprime(p2)}")
    print("-" * 20)

    # --- primes(n) demonstration ---
    n = 100
    prime_list = prime_checker.primes(n)
    print(f"Primes up to {n}: {prime_list}")
    print("-" * 20)

    # --- Maximal size calculation ---
    calculate_max_size()

<<<1250>>>