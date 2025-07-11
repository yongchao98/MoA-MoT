import sys
import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations, using
    a pre-computed bit array (sieve). It is designed to satisfy the problem's
    constraints on time complexity and memory.
    """
    def __init__(self, max_val=10000):
        self.max_val = max_val
        # The core data structure is a bit array, implemented with bytearray.
        # We need 1 bit for each number from 0 to max_val-1.
        # The size is ceil(max_val / 8).
        sieve_size = (max_val + 7) // 8
        self._sieve = bytearray([0xff] * sieve_size) # Initialize all bits to 1 (prime)

        # Mark 0 and 1 as not prime
        self._set_bit(0, False)
        self._set_bit(1, False)

        # Pre-compute primes using Sieve of Eratosthenes
        for p in range(2, int(math.sqrt(max_val))):
            if self.isprime(p):
                # Mark all multiples of p as not prime
                for i in range(p * p, max_val, p):
                    self._set_bit(i, False)

    def _set_bit(self, n, is_prime):
        """Helper to set or clear the nth bit in the bytearray."""
        byte_index = n // 8
        bit_index = n % 8
        if is_prime:
            # Set the bit to 1
            self._sieve[byte_index] |= (1 << bit_index)
        else:
            # Clear the bit to 0
            self._sieve[byte_index] &= ~(1 << bit_index)

    def isprime(self, p):
        """Checks if p is prime in O(1) time."""
        if not 0 <= p < self.max_val:
            return False # Out of bounds
        byte_index = p // 8
        bit_index = p % 8
        return (self._sieve[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """Returns a list of all primes <= n in O(n) time."""
        if not 0 <= n < self.max_val:
            raise ValueError(f"Input n must be between 0 and {self.max_val-1}")
        return [i for i in range(n + 1) if self.isprime(i)]

# --- Size Calculation ---
# The maximal size of the data structure is determined by the size of the
# bit array required to store the primality for all numbers up to the limit.

# The limit for p and n is < 10000, so we need to handle numbers 0-9999.
limit = 10000

# There are 8 bits in one byte.
bits_per_byte = 8

# The calculation for the required size in bytes is the ceiling of
# (total numbers / bits per byte).
# In integer arithmetic, ceiling(a / b) can be calculated as (a + b - 1) // b.
data_structure_size = (limit + bits_per_byte - 1) // bits_per_byte

print("--- Data Structure Size Calculation ---")
print("The design uses a bit array to store the primality for numbers up to 9999.")
print("\nThe final equation for the size in bytes is: ceil(total_numbers / bits_per_byte)")
print(f"Total numbers to store flags for = {limit}")
print(f"Number of bits in one byte = {bits_per_byte}")
print(f"Size = ({limit} + {bits_per_byte} - 1) // {bits_per_byte}")
print(f"The maximal size of this data structure is {data_structure_size} bytes.")