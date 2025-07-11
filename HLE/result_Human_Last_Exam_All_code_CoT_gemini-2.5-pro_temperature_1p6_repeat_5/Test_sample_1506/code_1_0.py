import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations using a
    pre-computed bit array populated by the Sieve of Eratosthenes.
    """

    def __init__(self, limit=10000):
        """
        Initializes the data structure by pre-computing primes up to the limit.
        """
        if limit < 2:
            raise ValueError("Limit must be at least 2.")
        self.limit = limit
        # Each byte holds 8 bits, so we need ceil(limit / 8) bytes.
        # The (limit + 7) // 8 is a trick for integer ceiling division.
        self._size_in_bytes = (limit + 7) // 8
        
        # Initialize bytearray with all bits set to 1 (all numbers are assumed prime initially).
        # 0xff is hexadecimal for 11111111 in binary.
        self._prime_bits = bytearray([0xff] * self._size_in_bytes)

        # --- Sieve of Eratosthenes Algorithm ---
        
        # 0 and 1 are not prime. Clear their corresponding bits.
        # Clear bit 0
        self._prime_bits[0 >> 3] &= ~(1 << (0 & 7))
        # Clear bit 1
        self._prime_bits[1 >> 3] &= ~(1 << (1 & 7))

        # Mark all multiples of primes as not prime (composite).
        for p in range(2, int(math.sqrt(limit)) + 1):
            # Check if p is still marked as prime.
            is_p_prime = (self._prime_bits[p >> 3] & (1 << (p & 7))) != 0
            if is_p_prime:
                # Mark all multiples of p starting from p*p.
                for i in range(p * p, limit, p):
                    # Clear the bit for number i.
                    self._prime_bits[i >> 3] &= ~(1 << (i & 7))

    def isprime(self, p):
        """
        Checks if p is a prime number in O(1) time.
        
        :param p: An integer where 0 <= p < 10000.
        :return: True if p is prime, False otherwise.
        """
        if not 0 <= p < self.limit:
            raise ValueError(f"Input p must be between 0 and {self.limit-1}")
        
        # Check the bit for number p. This is a constant time lookup.
        # p >> 3 is integer division p // 8 to find the byte index.
        # p & 7 is modulo p % 8 to find the bit index within the byte.
        return (self._prime_bits[p >> 3] & (1 << (p & 7))) != 0

    def primes(self, n):
        """
        Returns a list of all prime numbers p <= n in O(n) time.
        
        :param n: An integer where 0 <= n < 10000.
        :return: A list of prime numbers.
        """
        if not 0 <= n < self.limit:
            raise ValueError(f"Input n must be between 0 and {self.limit-1}")
        
        # Iterate from 2 to n and collect primes using the O(1) isprime check.
        return [i for i in range(2, n + 1) if self.isprime(i)]
    
    def get_max_size_in_bytes(self):
        """
        Returns the maximal size of the underlying bit array in bytes.
        """
        return self._size_in_bytes

# Main execution block
if __name__ == '__main__':
    # Create the data structure instance. The Sieve runs on initialization.
    prime_checker = PrimeDataStructure(limit=10000)

    # --- Maximal Size Calculation ---
    # The maximal size of the data structure is the size of its core component,
    # the bytearray, which is constant after initialization.
    
    total_numbers = prime_checker.limit
    bits_per_byte = 8
    max_size = prime_checker.get_max_size_in_bytes()

    print("To store the primality for each number from 0 to 9999, we need a bit array.")
    print(f"Total numbers to check: {total_numbers}")
    print("Each number's primality can be stored in a single bit.")
    print(f"The final equation for the size in bytes is: ceil(Total Numbers / Bits per Byte)")
    print(f"ceil({total_numbers} / {bits_per_byte}) = {max_size}")
    print(f"\nThe maximal size of this data structure is {max_size} bytes.")
