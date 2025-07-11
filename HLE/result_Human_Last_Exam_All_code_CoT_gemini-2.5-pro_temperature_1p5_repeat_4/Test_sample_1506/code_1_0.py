import math
import sys

class PrimeSieve:
    """
    A memory-efficient data structure for prime number operations using a bit array.
    """
    def __init__(self, max_num=9999):
        """
        Initializes the sieve up to max_num.
        """
        self.limit = max_num
        
        # Calculate the size needed for the bytearray.
        # We need one bit for each number from 0 to limit.
        num_of_bits = self.limit + 1
        self.size_in_bytes = (num_of_bits + 7) // 8
        
        # Initialize bytearray with all bits set to 1 (True), assuming all numbers are prime.
        self.sieve = bytearray([0xff] * self.size_in_bytes)

        # 0 and 1 are not prime; clear their bits.
        self._set_bit(0, 0)
        self._set_bit(1, 0)

        # Perform the Sieve of Eratosthenes
        for p in range(2, int(math.sqrt(self.limit)) + 1):
            # If p's bit is still set, it's a prime
            if self.isprime(p):
                # Mark all multiples of p as not prime
                # Start from p*p, as smaller multiples are already marked.
                for i in range(p * p, self.limit + 1, p):
                    self._set_bit(i, 0)

    def _set_bit(self, n, value):
        """Helper to set or clear a bit for number n."""
        byte_index = n >> 3  # n // 8
        bit_index = n & 7    # n % 8
        if value == 1:
            self.sieve[byte_index] |= (1 << bit_index)
        else:
            self.sieve[byte_index] &= ~(1 << bit_index)

    def isprime(self, p):
        """
        Checks if p is a prime number in O(1) time.
        p must be within the initialized limit.
        """
        if not 0 <= p <= self.limit:
            raise ValueError(f"Input {p} is out of the supported range [0, {self.limit}]")
        
        byte_index = p >> 3  # p // 8
        bit_index = p & 7    # p % 8
        
        # Return True if the bit is set (1), False otherwise (0).
        return (self.sieve[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all prime numbers p <= n in O(n) time.
        n must be within the initialized limit.
        """
        if not 0 <= n <= self.limit:
            raise ValueError(f"Input {n} is out of the supported range [0, {self.limit}]")

        # List comprehension iterates n+1 times, with an O(1) check inside.
        return [i for i in range(n + 1) if self.isprime(i)]

    def get_max_size_in_bytes(self):
        """Returns the size of the underlying data structure in bytes."""
        return self.size_in_bytes
        
    def print_size_calculation(self):
        """Prints the logic for the size calculation."""
        num_of_bits = self.limit + 1
        print("Size calculation:")
        print(f"Maximal number to check: {self.limit}")
        print(f"Total numbers to store (0 to {self.limit}): {num_of_bits}")
        print("Each number's primality is stored in a single bit.")
        print("There are 8 bits in 1 byte.")
        print("Final Equation for size in bytes: ceil(number_of_bits / bits_per_byte)")
        print(f"Size = ceil({num_of_bits} / 8) = {self.get_max_size_in_bytes()} bytes")


if __name__ == '__main__':
    # Create the data structure for primes up to 9999
    prime_checker = PrimeSieve(max_num=9999)

    # --- Demonstrate isprime(p) operation ---
    print("--- isprime(p) Demonstration ---")
    p1 = 29
    p2 = 30
    print(f"Is {p1} prime? {prime_checker.isprime(p1)}")
    print(f"Is {p2} prime? {prime_checker.isprime(p2)}")
    print(f"Is 9973 (largest prime < 10000) prime? {prime_checker.isprime(9973)}")
    print("-" * 20)

    # --- Demonstrate primes(n) operation ---
    print("--- primes(n) Demonstration ---")
    n = 30
    prime_list = prime_checker.primes(n)
    print(f"Primes up to {n}: {prime_list}")
    print("-" * 20)

    # --- Calculate and print the maximal size ---
    print("--- Data Structure Size ---")
    prime_checker.print_size_calculation()
    print("-" * 20)
    
    # Final answer for the maximal size
    print(f"The maximal size of this data structure is: {prime_checker.get_max_size_in_bytes()} bytes")
