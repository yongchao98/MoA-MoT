import sys
import math

class PrimeDataStructure:
    """
    An efficient data structure for prime number operations using a bit array.
    
    This data structure pre-computes prime numbers up to a specified limit using
    the Sieve of Eratosthenes and stores the results in a memory-efficient
    bit array (implemented with a bytearray).
    """

    def __init__(self, limit=10000):
        """
        Initializes the sieve for numbers up to the limit.
        """
        if limit <= 1:
            raise ValueError("Limit must be greater than 1.")
        self.limit = limit
        # Each bit represents a number, so we need limit bits.
        # A bytearray is used for memory efficiency.
        sieve_size_in_bytes = (limit + 7) // 8
        self._sieve = bytearray([0xFF] * sieve_size_in_bytes) # Initialize all to 1 (prime)

        # Helper to clear a bit (mark a number as not prime)
        def _clear_bit(n):
            byte_index = n >> 3  # n // 8
            bit_index = n & 7   # n % 8
            self._sieve[byte_index] &= ~(1 << bit_index)

        # Helper to get a bit's value
        def _get_bit(n):
            byte_index = n >> 3
            bit_index = n & 7
            return (self._sieve[byte_index] >> bit_index) & 1

        # 0 and 1 are not prime numbers
        _clear_bit(0)
        _clear_bit(1)

        # Sieve of Eratosthenes algorithm
        for number in range(2, int(math.sqrt(limit)) + 1):
            if _get_bit(number):  # If number is prime
                # Mark all its multiples as not prime
                for multiple in range(number * number, limit, number):
                    _clear_bit(multiple)

    def isprime(self, p):
        """
        Checks if p is a prime number in O(1) time.
        p must be less than the initialized limit.
        """
        if not 0 <= p < self.limit:
            raise ValueError(f"Input p must be within the range [0, {self.limit - 1}]")
        
        byte_index = p >> 3
        bit_index = p & 7
        return ((self._sieve[byte_index] >> bit_index) & 1) == 1

    def primes(self, n):
        """
        Returns a list of all primes p <= n in O(n) time.
        n must be less than the initialized limit.
        """
        if not 0 <= n < self.limit:
            raise ValueError(f"Input n must be within the range [0, {self.limit - 1}]")
        
        return [i for i in range(n + 1) if self.isprime(i)]

    def get_data_structure_size(self):
        """
        Returns the size of the core bytearray in bytes.
        """
        return len(self._sieve)

if __name__ == "__main__":
    # The problem specifies a limit of 10000
    LIMIT = 10000
    
    # Instantiate the data structure
    prime_checker = PrimeDataStructure(limit=LIMIT)

    print("--- Data Structure Functionality Demo ---")
    print(f"Is 9973 prime? {prime_checker.isprime(9973)}")
    print(f"Is 100 prime? {prime_checker.isprime(100)}")
    print(f"Primes up to 30: {prime_checker.primes(30)}\n")

    # --- Maximal Size Calculation ---
    print("--- Maximal Size Calculation ---")
    
    # We need to store primality for numbers from 0 to 9999.
    numbers_to_store = 10000
    bits_per_byte = 8
    
    # Calculate the size in bytes using integer arithmetic for ceiling division.
    # This matches how the bytearray size was calculated.
    size_in_bytes = (numbers_to_store + bits_per_byte - 1) // bits_per_byte
    
    print(f"The data structure needs to store primality for {numbers_to_store} numbers (0 to 9999).")
    print(f"Using a bit array, this requires {numbers_to_store} bits.")
    print("The size in bytes is calculated by dividing the number of bits by the bits in a byte (8).")
    print("Final Equation:")
    print(f"ceil({numbers_to_store} / {bits_per_byte}) = {size_in_bytes} bytes")
    
    # Verify with the actual created object's size
    actual_size = prime_checker.get_data_structure_size()
    print(f"\nThe actual size of the implemented data structure is {actual_size} bytes.")

<<<1250>>>