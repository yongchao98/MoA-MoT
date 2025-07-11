import math

class PrimeChecker:
    """
    A memory-efficient data structure for prime number operations up to a limit,
    using a Sieve of Eratosthenes implemented with a bit array (bytearray).
    """

    def __init__(self, limit=10000):
        """
        Initializes the Sieve of Eratosthenes up to the given limit.
        This pre-computation step populates the sieve.
        """
        if limit < 2:
            raise ValueError("Limit must be at least 2.")
        self.limit = limit
        # Calculate the size of the bytearray needed. (limit + 7) // 8 is a
        # common way to write ceiling division for integers.
        sieve_size = (limit + 7) // 8
        # Initialize a bytearray with all bits set to 1 (0xFF), marking all
        # numbers as potentially prime.
        self.sieve = bytearray([0xFF] * sieve_size)

        # Helper functions for bit manipulation for code clarity.
        # These will be used to implement the Sieve algorithm.
        def _clear_bit(k):
            # Sets the bit corresponding to number k to 0 (not prime).
            byte_index = k // 8
            bit_index = k % 8
            # The mask `~(1 << bit_index)` has a 0 at the target bit position
            # and 1s elsewhere. ANDing it with the byte clears the bit.
            self.sieve[byte_index] &= ~(1 << bit_index)

        def _is_set(k):
            # Checks if the bit corresponding to number k is 1 (prime).
            byte_index = k // 8
            bit_index = k % 8
            # Shift the bit to the rightmost position and check if it's 1.
            return (self.sieve[byte_index] >> bit_index) & 1

        # Mark 0 and 1 as not prime.
        _clear_bit(0)
        _clear_bit(1)

        # Sieve of Eratosthenes algorithm.
        # We only need to check for primes up to sqrt(limit).
        for p in range(2, int(math.sqrt(limit)) + 1):
            if _is_set(p):  # If p is prime...
                # ...then mark all of its multiples as not prime.
                # We can start from p*p, as smaller multiples (like 2*p, 3*p)
                # would have been marked by smaller primes (2, 3, etc.).
                for i in range(p * p, limit, p):
                    _clear_bit(i)

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be within the range [0, limit-1].
        """
        if not 0 <= p < self.limit:
            raise ValueError(f"Input p must be between 0 and {self.limit - 1}")
        
        byte_index = p // 8
        bit_index = p % 8
        return (self.sieve[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all prime numbers p <= n in O(n) time.
        n must be within the range [0, limit-1].
        """
        if not 0 <= n < self.limit:
            raise ValueError(f"Input n must be between 0 and {self.limit - 1}")
        
        # Iterate from 2 to n and use the O(1) isprime check.
        return [i for i in range(2, n + 1) if self.isprime(i)]

    def get_data_structure_size_in_bytes(self):
        """
        Returns the size of the core data structure (the bytearray).
        """
        return len(self.sieve)

# --- Main execution ---
LIMIT = 10000
checker = PrimeChecker(limit=LIMIT)

# Calculate and print the maximal size of the data structure.
NUM_OF_ITEMS = LIMIT
BITS_PER_BYTE = 8
size_in_bytes = checker.get_data_structure_size_in_bytes()

print("Calculating the maximal size of the data structure:")
print("The data structure needs to store the primality for each number from 0 to 9999.")
# As requested, outputting the numbers in the final equation.
print(f"Total numbers to check: {NUM_OF_ITEMS}")
print(f"Bits per number: 1")
print(f"Bits per byte: {BITS_PER_BYTE}")
print(f"Equation for size in bytes: ceil({NUM_OF_ITEMS} / {BITS_PER_BYTE})")
print(f"Resulting size = {size_in_bytes} bytes")
print(f"\nThe maximal size of this data structure in bytes is: {size_in_bytes}")

<<<1250>>>