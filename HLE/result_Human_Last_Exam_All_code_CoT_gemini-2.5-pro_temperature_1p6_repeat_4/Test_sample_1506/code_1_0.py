import math

class PrimeSieve:
    """
    A memory-efficient data structure for prime number operations using a
    pre-computed Sieve of Eratosthenes stored in a bit array.

    - isprime(p): O(1) - Checks if p is a prime.
    - primes(n):  O(n) - Returns a list of primes up to n.
    """

    def __init__(self, max_num=9999):
        """
        Initializes the sieve for all numbers up to max_num.
        """
        if not isinstance(max_num, int) or max_num < 0:
            raise ValueError("max_num must be a non-negative integer.")
            
        self._max_num = max_num
        # Determine the size of the bytearray needed. We need max_num + 1 bits.
        # This is a robust way to calculate ceil((max_num + 1) / 8).
        self._sieve_size_in_bytes = (self._max_num + 8) // 8
        
        # Create a bytearray, initializing all bits to 1 (True).
        self._sieve = bytearray([0xFF] * self._sieve_size_in_bytes)

        # Helper to clear a bit (set to 0 for non-prime).
        def _clear_bit(n):
            byte_index = n >> 3  # equivalent to n // 8
            bit_index = n & 7    # equivalent to n % 8
            # Use a mask to clear the specific bit without affecting others.
            self._sieve[byte_index] &= ~(1 << bit_index)

        # 0 and 1 are not prime numbers.
        if self._max_num >= 0: _clear_bit(0)
        if self._max_num >= 1: _clear_bit(1)
        
        # Populate the sieve by marking non-primes.
        # We only need to iterate up to the square root of the maximum number.
        for number in range(2, int(math.sqrt(self._max_num)) + 1):
            # Check if the number is prime (its bit is still 1).
            if (self._sieve[number >> 3] >> (number & 7)) & 1:
                # Mark all multiples of the number as not prime.
                # Start marking from number*number.
                for multiple in range(number * number, self._max_num + 1, number):
                    _clear_bit(multiple)

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be within the pre-computed range [0, max_num].
        """
        if not (0 <= p <= self._max_num):
            raise ValueError(f"Input {p} is out of the supported range [0, {self._max_num}]")
        
        # Find the byte and the bit for the number p.
        byte_index = p >> 3
        bit_index = p & 7
        
        # Check if the bit at the position is set to 1.
        return (self._sieve[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all prime numbers <= n in O(n) time.
        n must be within the pre-computed range [0, max_num].
        """
        if not (0 <= n <= self._max_num):
            raise ValueError(f"Input {n} is out of the supported range [0, {self._max_num}]")
        
        # Iterate from 0 to n and collect numbers that are prime.
        return [i for i in range(n + 1) if self.isprime(i)]

    def get_data_structure_size(self):
        """Returns the size of the core data structure in bytes."""
        return self._sieve_size_in_bytes

if __name__ == '__main__':
    # The problem specifies that p and n will be less than 10000.
    # Therefore, we need to handle numbers up to 9999.
    MAX_NUMBER = 9999

    # The number of boolean flags we need to store (for numbers 0 to 9999).
    numbers_to_store = MAX_NUMBER + 1
    bits_per_byte = 8

    # Calculate the size in bytes.
    # We use integer division and add 1 if there's a remainder, or use ceiling division.
    # (numbers_to_store + bits_per_byte - 1) // bits_per_byte is a common way to do ceiling division.
    size_in_bytes = (numbers_to_store + bits_per_byte - 1) // bits_per_byte
    
    print("This program designs a memory-efficient data structure for prime lookups.")
    print("\n--- Calculation of Maximal Data Structure Size ---")
    print(f"To handle numbers less than 10000, we need to store primality for numbers from 0 to {MAX_NUMBER}.")
    print(f"Total boolean flags required: {numbers_to_store}")
    print(f"Using a bit array, each flag takes 1 bit.")
    print(f"There are {bits_per_byte} bits in 1 byte.")
    print("\nFinal equation for the size in bytes:")
    print(f"size = ceil({numbers_to_store} / {bits_per_byte})")
    print(f"size = {size_in_bytes} bytes")

    # Example Usage (optional, uncomment to test)
    # sieve = PrimeSieve(MAX_NUMBER)
    # print("\n--- Example Usage ---")
    # print(f"Is 9973 (a prime) prime? {sieve.isprime(9973)}")
    # print(f"Is 100 (not a prime) prime? {sieve.isprime(100)}")
    # print(f"Primes up to 20: {sieve.primes(20)}")

<<<1250>>>