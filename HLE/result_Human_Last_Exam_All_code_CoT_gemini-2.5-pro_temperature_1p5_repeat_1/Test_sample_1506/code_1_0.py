import math
import sys

class PrimeChecker:
    """
    A memory-efficient data structure for prime number operations.
    It uses a pre-computed bit array (sieve) to achieve fast lookups.
    """
    def __init__(self, max_num=9999):
        """
        Initializes the PrimeChecker by pre-computing primes up to max_num
        using a Sieve of Eratosthenes and storing them in a bit array.
        """
        if max_num >= 10000:
            raise ValueError("This class is designed for numbers less than 10000.")
            
        self.max_num = max_num
        # A bit array to store primality, implemented with a bytearray.
        # We need to store primality for numbers 0 to 9999 (10000 total).
        # Size in bytes = ceil(10000 / 8) = 1250.
        sieve_size = (max_num + 8) // 8
        self._sieve = bytearray([0b11111111] * sieve_size)

        # Helper to set a bit to 0 (mark as not prime)
        def _clear_bit(n):
            byte_index = n // 8
            bit_index = n % 8
            self._sieve[byte_index] &= ~(1 << bit_index)

        # 0 and 1 are not prime
        _clear_bit(0)
        _clear_bit(1)

        # Sieve of Eratosthenes algorithm
        for p in range(2, int(math.sqrt(max_num)) + 1):
            # Check if p is prime (its bit is still 1)
            if self.isprime(p):
                # Mark all multiples of p (starting from p*p) as not prime
                for i in range(p * p, max_num + 1, p):
                    _clear_bit(i)

    def isprime(self, p: int) -> bool:
        """
        Checks if p is a prime number in O(1) time.
        p must be less than 10000.
        """
        if not (0 <= p <= self.max_num):
            return False 
        
        byte_index = p // 8
        bit_index = p % 8
        # O(1) lookup and bitwise operation
        return (self._sieve[byte_index] >> bit_index) & 1 == 1

    def primes(self, n: int) -> list:
        """
        Returns a list of all primes p <= n in O(n) time.
        n must be less than 10000.
        """
        if not (0 <= n <= self.max_num):
             raise ValueError(f"Input n must be between 0 and {self.max_num}.")
        
        # Iterate from 0 to n (O(n)) and perform an O(1) check for each.
        return [i for i in range(n + 1) if self.isprime(i)]

    def get_data_structure_size_in_bytes(self) -> int:
        """
        Returns the size of the core data structure (_sieve) in bytes.
        """
        return len(self._sieve)

# Main execution block to demonstrate usage and calculate size
if __name__ == "__main__":
    # Create an instance of the PrimeChecker for numbers up to 9999.
    prime_checker = PrimeChecker(max_num=9999)

    print("--- Demonstrating Operations ---")
    
    # 1. isprime(p): O(1)
    p1 = 9973 # A prime number
    p2 = 9974 # Not a prime number
    print(f"Is {p1} prime? {prime_checker.isprime(p1)}")
    print(f"Is {p2} prime? {prime_checker.isprime(p2)}")
    
    print("-" * 20)

    # 2. primes(n): O(n)
    n = 30
    prime_list = prime_checker.primes(n)
    print(f"Primes less than or equal to {n}: {prime_list}")

    print("-" * 20)
    
    # --- Calculate and print the data structure size ---
    print("--- Calculating Data Structure Size ---")
    
    total_numbers = prime_checker.max_num + 1
    bits_per_byte = 8
    size_in_bytes = prime_checker.get_data_structure_size_in_bytes()

    print(f"Total numbers covered: {total_numbers} (from 0 to 9999)")
    print(f"Bits used per number: 1")
    print(f"Total bits required: {total_numbers}")
    print(f"Calculation: {total_numbers} bits / {bits_per_byte} bits/byte = {size_in_bytes} bytes")
    print(f"Maximal size of this data structure in bytes: {size_in_bytes}")
