import math

class PrimeDataStructure:
    """
    A data structure for efficient prime number operations using a pre-computed sieve.
    The sieve is generated once upon initialization.
    """
    def __init__(self, limit=10000):
        """
        Initializes the data structure by generating a prime sieve up to the limit
        using the Sieve of Eratosthenes algorithm.
        """
        self.limit = limit
        # For implementation clarity, we use a list of booleans.
        # A true bit array would be the most memory-efficient.
        self._sieve = [True] * limit
        if limit > 1:
            self._sieve[0] = self._sieve[1] = False
        
        # Mark all multiples of primes as not prime
        for i in range(2, int(math.sqrt(limit)) + 1):
            if self._sieve[i]:
                for multiple in range(i * i, limit, i):
                    self._sieve[multiple] = False

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time by looking it up in the sieve.
        p must be less than the initialized limit.
        """
        if 0 <= p < self.limit:
            return self._sieve[p]
        # Handle out-of-bounds input
        raise ValueError(f"Input {p} is out of the supported range [0, {self.limit-1}]")

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n in O(n) time.
        n must be less than the initialized limit.
        """
        if not 0 <= n < self.limit:
             raise ValueError(f"Input {n} is out of the supported range [0, {self.limit-1}]")
            
        # Iterate from 2 to n and collect primes
        return [i for i in range(2, n + 1) if self._sieve[i]]

# --- Main execution ---

# 1. Create an instance of the data structure.
prime_checker = PrimeDataStructure(limit=10000)

# 2. Demonstrate the 'isprime' operation (O(1)).
p1 = 9973 # A known prime
p2 = 9999 # Not a prime
print(f"Is {p1} prime? {prime_checker.isprime(p1)}")
print(f"Is {p2} prime? {prime_checker.isprime(p2)}")

# 3. Demonstrate the 'primes' operation (O(n)).
n = 30
print(f"Primes up to {n}: {prime_checker.primes(n)}")

# 4. Calculate and print the maximal size of the data structure in bytes.
# The most memory-efficient implementation is a bit array.
# We need to store primality for 10,000 numbers (0 to 9999).
# This requires 10,000 bits. There are 8 bits in a byte.
num_bits = 10000
bits_in_byte = 8
max_size_bytes = num_bits / bits_in_byte

print("\n--- Maximal Size Calculation ---")
print("The most efficient data structure is a bit array.")
print(f"It needs to store primality for {num_bits} numbers, requiring {num_bits} bits.")
print(f"The size in bytes is the number of bits divided by the number of bits in a byte ({bits_in_byte}).")
print(f"Calculation: {num_bits} / {bits_in_byte} = {int(max_size_bytes)}")
print(f"Maximal size of the data structure: {int(max_size_bytes)} bytes.")

<<<1250>>>