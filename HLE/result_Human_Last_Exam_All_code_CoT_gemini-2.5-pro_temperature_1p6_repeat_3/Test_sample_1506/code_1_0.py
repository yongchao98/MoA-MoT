import sys
import math

# Plan:
# 1. Design a memory-efficient data structure for prime number lookups.
#    - To achieve O(1) for `isprime`, a pre-computed lookup table is ideal.
#    - To be memory-efficient, a bit array is the best choice. Python's `bytearray` can implement this.
#    - The "Sieve of Eratosthenes" algorithm is perfect for pre-computing the primes and populating the bit array.
#
# 2. Implement the data structure in a Python class (`PrimeDataStructure`).
#    - The constructor `__init__` will build the sieve for numbers up to the required limit (9999).
#    - `isprime(p)` will perform a bitwise check on the `bytearray` (O(1)).
#    - `primes(n)` will iterate from 2 to `n`, using `isprime()` for each number (O(n)).
#
# 3. Calculate the maximal size of the data structure.
#    - The dominant part of the structure is the `bytearray`.
#    - The size of the data is `ceil((limit + 1) / 8)` bytes.
#    - The total size in memory, including Python's object overhead, is found using `sys.getsizeof()`.
#    - The script will calculate and print this size with a clear breakdown.

class PrimeDataStructure:
    """
    An efficient data structure for prime number operations using a bit-array sieve.
    """

    def __init__(self, limit: int):
        """
        Initializes the sieve for all numbers up to the specified limit.
        """
        if limit < 2:
            raise ValueError("Limit must be at least 2.")
        self.limit = limit
        
        # Calculate size for the bytearray. Each byte stores 8 numbers' primality.
        # We need to store primality for numbers from 0 to `limit`.
        self.sieve_size_in_bytes = (limit + 1 + 7) // 8
        
        # Initialize bytearray. All bits are set to 1 (True, i.e., potentially prime).
        self._sieve = bytearray(b'\xff' * self.sieve_size_in_bytes)

        # Helper to set a bit to 0 (False, i.e., not prime)
        def _clear_bit(n):
            byte_index = n >> 3  # n // 8
            bit_index = n & 7   # n % 8
            self._sieve[byte_index] &= ~(1 << bit_index)

        # 0 and 1 are not prime numbers.
        _clear_bit(0)
        _clear_bit(1)

        # Sieve of Eratosthenes algorithm
        for p in range(2, int(math.sqrt(limit)) + 1):
            if self.isprime(p):
                for i in range(p * p, limit + 1, p):
                    _clear_bit(i)

    def isprime(self, p: int) -> bool:
        """
        Checks if a number p is prime in O(1) time.
        """
        if not 0 <= p <= self.limit:
            raise ValueError(f"Input must be between 0 and {self.limit}.")
        
        byte_index = p >> 3
        bit_index = p & 7
        return (self._sieve[byte_index] >> bit_index) & 1 == 1

    def primes(self, n: int) -> list[int]:
        """
        Returns a list of all primes <= n in O(n) time.
        """
        if not 0 <= n <= self.limit:
            raise ValueError(f"Input must be between 0 and {self.limit}.")
        
        return [i for i in range(n + 1) if self.isprime(i)]

    def get_data_structure_size(self) -> int:
        """
        Returns the size of the core data structure (_sieve) in bytes.
        """
        return sys.getsizeof(self._sieve)

# Set the limit based on the problem statement (p < 10000, n < 10000)
# This means we need to handle numbers up to 9999.
LIMIT = 9999

# Create an instance of the data structure.
prime_checker = PrimeDataStructure(limit=LIMIT)

# Demonstrate the 'isprime' operation
p1 = 9973  # A known prime
print(f"Is {p1} a prime? {prime_checker.isprime(p1)}")
p2 = 9998  # A non-prime
print(f"Is {p2} a prime? {prime_checker.isprime(p2)}\n")

# Demonstrate the 'primes' operation
n_limit = 30
print(f"Primes up to {n_limit}: {prime_checker.primes(n_limit)}\n")

# --- Calculate and print the maximal size of the data structure ---

print("--- Data Structure Size Calculation ---")

# The number of integers for which we need to store primality (0 to 9999)
num_of_integers = LIMIT + 1

# Calculate the number of bytes needed for the data payload
data_size_bytes = prime_checker.sieve_size_in_bytes

# Get the total size of the bytearray object in memory
total_size_bytes = prime_checker.get_data_structure_size()

# The overhead is the difference between the total size and the data size
overhead_bytes = total_size_bytes - data_size_bytes

print(f"To support numbers < 10000, we must store primality for integers from 0 to {LIMIT}.")
print(f"Total integers to track: {num_of_integers}")
print(f"Size of data payload = ceil({num_of_integers} / 8) = {data_size_bytes} bytes")
print(f"Python object overhead = {overhead_bytes} bytes")
print("Maximal size of data structure = Overhead + Data Size")
print(f"Final Equation: {overhead_bytes} + {data_size_bytes} = {total_size_bytes}")

<<<1311>>>