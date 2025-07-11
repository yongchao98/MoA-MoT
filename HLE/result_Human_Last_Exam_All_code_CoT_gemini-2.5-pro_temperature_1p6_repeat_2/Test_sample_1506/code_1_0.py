import math
import sys

class PrimeChecker:
    """
    A memory-efficient data structure for prime number operations
    using a bit array implemented with a bytearray (Sieve of Eratosthenes).
    """

    def __init__(self, max_num):
        """
        Initializes the sieve for all numbers up to max_num.
        """
        if max_num < 0:
            raise ValueError("Maximum number must be non-negative.")
        
        self.max_num = max_num
        # Size of the bit array in bytes
        sieve_size_bytes = (max_num + 7) // 8
        
        # The core data structure: a bytearray acting as a bit array.
        # Initialize all bits to 1 (True).
        self.sieve = bytearray([0b11111111]) * sieve_size_bytes
        
        # Mark 0 and 1 as not prime by clearing their corresponding bits.
        self._clear_bit(0)
        if max_num >= 1:
            self._clear_bit(1)
        
        # Perform the Sieve of Eratosthenes
        for p in range(2, int(math.sqrt(max_num)) + 1):
            if self._is_bit_set(p): # If p is prime
                # Mark all multiples of p as not prime
                for i in range(p * p, max_num, p):
                    self._clear_bit(i)

    def _is_bit_set(self, n):
        """Checks if the bit for number n is set (is prime)."""
        byte_index = n // 8
        bit_index = n % 8
        return (self.sieve[byte_index] & (1 << bit_index)) != 0

    def _clear_bit(self, n):
        """Clears the bit for number n (marks as not prime)."""
        byte_index = n // 8
        bit_index = n % 8
        self.sieve[byte_index] &= ~(1 << bit_index)

    def isprime(self, p):
        """
        Checks if a number p is prime.
        Time Complexity: O(1).
        """
        if not 0 <= p < self.max_num:
            raise ValueError(f"Input {p} is out of the pre-computed range [0, {self.max_num-1}]")
        return self._is_bit_set(p)

    def primes(self, n):
        """
        Returns a list of all prime numbers less than or equal to n.
        Time Complexity: O(n).
        """
        if not 0 <= n < self.max_num:
            raise ValueError(f"Input {n} is out of the pre-computed range [0, {self.max_num-1}]")
        
        prime_list = []
        for i in range(2, n + 1):
            if self.isprime(i):
                prime_list.append(i)
        return prime_list

if __name__ == "__main__":
    LIMIT = 10000
    prime_checker = PrimeChecker(LIMIT)

    print(f"Data structure designed: Bit Array (Sieve) up to {LIMIT-1}")
    print("Time complexity for isprime(p): O(1)")
    print("Time complexity for primes(n): O(n)\n")
    
    # --- Example Usage ---
    p1 = 9973 # A prime number
    p2 = 9974 # Not a prime number
    print(f"Is {p1} a prime number? {prime_checker.isprime(p1)}")
    print(f"Is {p2} a prime number? {prime_checker.isprime(p2)}\n")
    
    n = 100
    print(f"Primes up to {n}:")
    print(prime_checker.primes(n))
    print("\n" + "="*40)
    
    # --- Size Calculation ---
    data_structure_size_bytes = sys.getsizeof(prime_checker.sieve)
    
    # Calculation explanation
    bytes_for_bits = (LIMIT + 7) // 8
    overhead = data_structure_size_bytes - bytes_for_bits
    
    print("Maximal Size Calculation:")
    print(f"The data structure is a bytearray representing {LIMIT} bits.")
    print(f"Number of bytes needed for bits = ceil({LIMIT} / 8) = {bytes_for_bits} bytes.")
    print(f"Python's bytearray object adds an overhead of {overhead} bytes.")
    print(f"Total Size = {bytes_for_bits} + {overhead} = {data_structure_size_bytes} bytes.")
    print("="*40)
    print(f"\nThe maximal size of this data structure is {data_structure_size_bytes} bytes.")
    
    # Final answer in the required format
    print(f"<<<{data_structure_size_bytes}>>>")