import math
import sys

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations
    using a pre-computed bit array via the Sieve of Eratosthenes.
    It is initialized for all numbers less than `max_num`.
    """
    def __init__(self, max_num=10000):
        """Initializes the bit array for numbers 0 to max_num-1."""
        if max_num < 0:
            raise ValueError("max_num must be non-negative")
        self.max_num = max_num
        
        # The core data structure: a bit array implemented with bytearray.
        # Size in bytes is ceiling(max_num / 8).
        self.size_in_bytes = (self.max_num + 7) // 8
        self._sieve_bits = bytearray([0xff] * self.size_in_bytes) # Assume all are prime (bit=1)

        # Helper functions for bit manipulation, defined locally for encapsulation
        def _get_bit(n):
            # Internal helper to get the primality bit for number n.
            byte_index = n >> 3 # n // 8
            bit_index = n & 7  # n % 8
            return (self._sieve_bits[byte_index] >> bit_index) & 1

        def _clear_bit(n):
            # Internal helper to set the primality bit for number n to 0 (not prime).
            byte_index = n >> 3 # n // 8
            bit_index = n & 7  # n % 8
            self._sieve_bits[byte_index] &= ~(1 << bit_index)

        # --- Sieve of Eratosthenes Initialization ---
        # 0 and 1 are not prime.
        if self.max_num > 0:
            _clear_bit(0)
        if self.max_num > 1:
            _clear_bit(1)

        # Iterate from 2 up to sqrt(max_num)
        for p in range(2, int(math.sqrt(self.max_num)) + 1):
            if _get_bit(p): # If p is prime
                # Mark all multiples of p as not prime.
                # Start from p*p, as smaller multiples have already been marked by smaller primes.
                for i in range(p * p, self.max_num, p):
                    _clear_bit(i)

    def isprime(self, p):
        """
        Checks if p is a prime number in O(1) time.
        p must be less than the initial max_num.
        """
        if not 0 <= p < self.max_num:
            raise ValueError(f"Input p must be between 0 and {self.max_num-1}")

        byte_index = p >> 3 # p // 8
        bit_index = p & 7  # p % 8
        return (self._sieve_bits[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all primes p <= n in O(n) time.
        n must be less than the initial max_num.
        """
        if not 0 <= n < self.max_num:
            raise ValueError(f"Input n must be between 0 and {self.max_num-1}")

        # The check for each number is O(1), and the loop runs n+1 times. Total time is O(n).
        return [i for i in range(n + 1) if self.isprime(i)]

# --- Main execution to calculate and print the size ---

# The problem specifies that the functions will handle inputs p < 10000 and n < 10000.
# This means our data structure must store primality information for numbers 0 through 9999.
UPPER_LIMIT = 10000

# The data is stored in a bit array, requiring 1 bit for each number.
total_bits = UPPER_LIMIT

# There are 8 bits in one byte.
bits_per_byte = 8

# To find the size in bytes, we divide the total bits by the number of bits in a byte.
# We must take the ceiling of the result, as we cannot allocate a fraction of a byte.
# The ceiling of (a / b) for positive integers can be calculated using integer division: (a + b - 1) // b.
size_in_bytes = (total_bits + bits_per_byte - 1) // bits_per_byte

print("Calculating the maximal size of the data structure in bytes:")
print(f"The structure must hold information for numbers from 0 to {UPPER_LIMIT - 1}.")
print(f"Total numbers to represent = {total_bits}")
print(f"Bits needed (1 per number) = {total_bits}")
print(f"Number of bits per byte = {bits_per_byte}")
print(f"Maximal size = ceil(Total bits / Bits per byte) = ceil({total_bits} / {bits_per_byte}) = {size_in_bytes} bytes")

# Note: This represents the fundamental size of the data itself. A Python object
# holding this data, such as the `bytearray` used in the class, will have some
# additional memory overhead for the object's metadata. The size of the core data
# is the most accurate answer to the design question.
# For example, `sys.getsizeof(bytearray(1250))` on a 64-bit system is typically 1307 bytes.
# However, the required information storage is 1250 bytes.

# Example usage (optional):
# prime_checker = PrimeDataStructure()
# print(f"Is 9973 a prime? {prime_checker.isprime(9973)}")
# print(f"Primes up to 30: {prime_checker.primes(30)}")