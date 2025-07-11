import math
import sys

class PrimeHandler:
    """
    A data structure for efficient prime number operations using a pre-computed bit sieve.
    This design is the most efficient in terms of memory.
    """
    def __init__(self, limit=9999):
        """
        Initializes the data structure by building a sieve up to the specified limit.
        """
        self.limit = limit
        # A bit array is the most memory-efficient structure. We use a Python bytearray.
        # The size is ceiling of (limit + 1) / 8.
        self.sieve_size_in_bytes = (limit + 8) // 8
        self._sieve = bytearray(b'\x01') * self.sieve_size_in_bytes # Initialize all bits to 1 (prime)

        # Helper functions to get/set bits in the bytearray
        self._is_bit_set = lambda n: (self._sieve[n >> 3] >> (n & 7)) & 1
        self._clear_bit = lambda n: self._sieve.__setitem__(n >> 3, self._sieve[n >> 3] & ~(1 << (n & 7)))

        # Sieve of Eratosthenes algorithm to pre-populate the bit array
        # Mark 0 and 1 as not prime
        self._clear_bit(0)
        self._clear_bit(1)

        # Mark all multiples of primes as not prime
        for num in range(2, int(math.sqrt(limit)) + 1):
            if self._is_bit_set(num):
                for multiple in range(num * num, limit + 1, num):
                    self._clear_bit(multiple)

    def isprime(self, p):
        """
        Checks if p is a prime number. Time complexity: O(1).
        p must be < 10000.
        """
        if not 0 <= p <= self.limit:
            raise ValueError(f"Input p must be between 0 and {self.limit}.")
        return self._is_bit_set(p) == 1

    def primes(self, n):
        """
        Returns a list of all prime numbers p <= n. Time complexity: O(n).
        n must be < 10000.
        """
        if not 0 <= n <= self.limit:
            raise ValueError(f"Input n must be between 0 and {self.limit}.")
        # The loop runs n+1 times, and each call to isprime is O(1).
        return [i for i in range(n + 1) if self.isprime(i)]

    def get_data_structure_size_bytes(self):
        """
        Returns the size in bytes of the underlying data structure (the sieve).
        """
        return len(self._sieve)

# --- Main execution block ---
if __name__ == "__main__":
    # The problem specifies p < 10000 and n < 10000, so the max number is 9999.
    MAX_LIMIT = 9999

    # Instantiate the data structure
    prime_handler = PrimeHandler(limit=MAX_LIMIT)

    # --- Demonstrate the operations ---
    print("--- Demonstrating Operations ---")
    print(f"isprime(9973): {prime_handler.isprime(9973)}")
    print(f"isprime(100): {prime_handler.isprime(100)}")
    print(f"primes(30): {prime_handler.primes(30)}")
    print("-" * 30)
    print()

    # --- Calculate and explain the maximal size ---
    print("--- Maximal Size Calculation ---")
    print("The most efficient data structure is a bit array where each bit represents a number's primality.")

    num_of_integers = MAX_LIMIT + 1
    bits_per_byte = 8
    size_in_bytes = prime_handler.get_data_structure_size_bytes()

    print(f"Numbers to check: 0 to {MAX_LIMIT} (a total of {num_of_integers} numbers).")
    print(f"Storage per number: 1 bit.")
    print(f"Total bits required: {num_of_integers} bits.")
    print(f"There are {bits_per_byte} bits in 1 byte.")

    print("\nThe final equation for the size in bytes is:")
    print(f"Maximal Size = ceil(Total Numbers / Bits per Byte)")
    print(f"Maximal Size = ceil({num_of_integers} / {bits_per_byte}) = {size_in_bytes} bytes")

<<<1250>>>