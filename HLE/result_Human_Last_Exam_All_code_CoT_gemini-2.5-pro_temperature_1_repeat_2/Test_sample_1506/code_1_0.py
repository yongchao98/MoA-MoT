import sys
import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations based on a
    pre-computed Sieve of Eratosthenes stored in a bit array.

    - isprime(p): O(1)
    - primes(n): O(n)
    """

    def __init__(self, max_limit=10000):
        """
        Initializes the data structure by generating a sieve of primes.
        """
        if max_limit < 2:
            raise ValueError("max_limit must be at least 2")
        self.max_limit = max_limit
        # Each number needs one bit. We use a bytearray for memory efficiency.
        # The size is ceiling(max_limit / 8).
        byte_size = (self.max_limit + 7) // 8
        # Initialize all bits to 1 (True), assuming all numbers are prime.
        self.sieve = bytearray([0b11111111]) * byte_size

        # Helper to set a bit for number n to 0 (False) or 1 (True).
        def _set_bit(n, value):
            byte_index = n // 8
            bit_index = n % 8
            if not value:
                self.sieve[byte_index] &= ~(1 << bit_index)
            else:
                 self.sieve[byte_index] |= (1 << bit_index)

        # Mark 0 and 1 as not prime
        _set_bit(0, False)
        if self.max_limit > 0:
            _set_bit(1, False)

        # Sieve of Eratosthenes algorithm
        for p in range(2, int(math.sqrt(self.max_limit)) + 1):
            if self.isprime(p):
                # Mark all multiples of p as not prime
                for i in range(p * p, self.max_limit, p):
                    _set_bit(i, False)

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be less than the pre-computed max_limit.
        """
        if not 0 <= p < self.max_limit:
            raise ValueError(f"Input p must be between 0 and {self.max_limit-1}")

        byte_index = p // 8
        bit_index = p % 8
        return (self.sieve[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n in O(n) time.
        n must be less than the pre-computed max_limit.
        """
        if not 0 <= n < self.max_limit:
             raise ValueError(f"Input n must be between 0 and {self.max_limit-1}")

        return [i for i in range(2, n + 1) if self.isprime(i)]

def calculate_and_print_size():
    """
    Calculates and prints the size of the designed data structure.
    """
    max_num_to_check = 9999
    num_of_entries = max_num_to_check + 1 # For numbers 0 to 9999
    bits_per_byte = 8

    # The size in bytes is the ceiling of the total bits divided by 8.
    size_in_bytes = (num_of_entries + bits_per_byte - 1) // bits_per_byte

    print("The most memory-efficient data structure is a bit array that pre-computes")
    print("primality for all numbers up to the required limit (9999).")
    print("\nCalculating the maximal size of this data structure:")
    print(f"- Total numbers to store primality for (0 to {max_num_to_check}): {num_of_entries}")
    print("- Bits required per number: 1")
    print(f"- Bytes required ({bits_per_byte} bits per byte): ceil({num_of_entries} / {bits_per_byte})")
    print("\nFinal Equation:")
    print(f"{num_of_entries} / {bits_per_byte} = {num_of_entries / bits_per_byte}")
    print(f"ceil({num_of_entries / bits_per_byte}) = {size_in_bytes}")
    print(f"\nThe maximal size of this data structure is {size_in_bytes} bytes.")


if __name__ == '__main__':
    # This part demonstrates the functionality and calculates the size.
    
    # 1. Create the data structure
    # This proves the design is functional.
    prime_ds = PrimeDataStructure(max_limit=10000)

    # 2. Calculate and print the size of the data structure.
    calculate_and_print_size()

    # You can uncomment the following lines to test the methods.
    # print("\n--- Testing the data structure ---")
    # print(f"Is 9973 prime? {prime_ds.isprime(9973)}") # 9973 is a prime
    # print(f"Is 9999 prime? {prime_ds.isprime(9999)}") # 9999 is not a prime
    # print(f"Primes up to 30: {prime_ds.primes(30)}")
<<<1250>>>