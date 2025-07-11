import math
import sys

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations up to a limit.
    
    It uses a bit array (implemented with Python's bytearray) to store a
    pre-computed Sieve of Eratosthenes. This design provides:
    - isprime(p): O(1) time complexity
    - primes(n): O(n) time complexity
    """
    def __init__(self, limit=10000):
        """
        Initializes the data structure by building the sieve.
        """
        if limit <= 1:
            raise ValueError("Limit must be greater than 1.")
        self.LIMIT = limit
        
        # Calculate the required size for the bit array in bytes.
        # Each number's primality is stored as a single bit.
        self._sieve_size_bytes = (self.LIMIT + 7) // 8
        
        # Initialize all bits to 1 (True, i.e., potentially prime).
        # A bytearray is a mutable sequence of integers in the range 0 <= x < 256.
        self._sieve = bytearray(b'\xff' * self._sieve_size_bytes)
        
        # 0 and 1 are not prime numbers, so clear their corresponding bits.
        self._clear_bit(0)
        self._clear_bit(1)
        
        # Implement the Sieve of Eratosthenes algorithm.
        # We only need to iterate up to the square root of the limit.
        for p in range(2, int(math.sqrt(self.LIMIT)) + 1):
            # If p is still marked as prime...
            if self._is_set(p):
                # ...then mark all of its multiples as not prime.
                # We can start marking from p*p, as smaller multiples
                # would have been marked by smaller primes.
                for i in range(p * p, self.LIMIT, p):
                    self._clear_bit(i)

    def _is_set(self, n):
        """Checks if the bit for number 'n' is set to 1 (is prime)."""
        byte_index = n >> 3  # Equivalent to n // 8
        bit_index = n & 7    # Equivalent to n % 8
        return (self._sieve[byte_index] & (1 << bit_index)) != 0

    def _clear_bit(self, n):
        """Clears the bit for number 'n' to 0 (is not prime)."""
        byte_index = n >> 3  # Equivalent to n // 8
        bit_index = n & 7    # Equivalent to n % 8
        # Use bitwise AND with the inverse of the bit mask to clear the bit.
        mask = ~(1 << bit_index)
        self._sieve[byte_index] &= mask

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be less than the initialized limit.
        """
        if not 0 <= p < self.LIMIT:
            raise ValueError(f"Input p must be between 0 and {self.LIMIT - 1}")
        return self._is_set(p)

    def primes(self, n):
        """
        Returns a list of all prime numbers p <= n in O(n) time.
        n must be less than the initialized limit.
        """
        if not 0 <= n < self.LIMIT:
            raise ValueError(f"Input n must be between 0 and {self.LIMIT - 1}")
        
        return [i for i in range(2, n + 1) if self.isprime(i)]
        
    def calculate_and_print_size(self):
        """Calculates and prints the size of the core data structure."""
        print("Calculating the maximal size of the data structure:")
        print("-" * 50)
        num_numbers = self.LIMIT
        bits_per_byte = 8
        
        # The core data is the bit array (sieve).
        # Size in bytes = total bits / bits per byte.
        data_size_bytes = self._sieve_size_bytes
        
        print(f"The structure must store primality for numbers up to {num_numbers - 1}.")
        print(f"Total numbers to store: {num_numbers}")
        print("We use a bit array, where 1 bit represents one number's primality.")
        print(f"Total bits required: {num_numbers}")
        print(f"Number of bits in a byte: {bits_per_byte}")
        print("\nThe size is calculated as follows:")
        print(f"Final Equation: {num_numbers} (bits) / {bits_per_byte} (bits per byte) = {data_size_bytes} bytes")
        print("-" * 50)
        print(f"The theoretical maximal size of this data structure is {data_size_bytes} bytes.")
        
        # Optional: Print the size of the Python bytearray object for context.
        actual_size = sys.getsizeof(self._sieve)
        print(f"\nNote: The actual size of the Python bytearray object in this environment is {actual_size} bytes, which includes implementation-specific overhead.")


# Create an instance of the data structure and calculate its size.
prime_checker = PrimeDataStructure(limit=10000)
prime_checker.calculate_and_print_size()

<<<1250>>>