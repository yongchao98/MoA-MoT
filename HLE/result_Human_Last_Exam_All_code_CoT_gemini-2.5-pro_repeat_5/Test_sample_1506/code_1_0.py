import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations using a bit sieve.
    """
    def __init__(self, max_val=9999):
        """
        Initializes the data structure by pre-computing primes up to max_val
        using the Sieve of Eratosthenes.
        """
        # The sieve needs to store primality for numbers from 0 to max_val.
        self._MAX_N = max_val + 1
        
        # Calculate the size of the bit array in bytes.
        # (self._MAX_N + 7) // 8 is an efficient way to calculate ceil(self._MAX_N / 8).
        sieve_size_bytes = (self._MAX_N + 7) // 8
        
        # Create a bytearray, initializing all bits to 1 (True), assuming all numbers are prime.
        # 0xff is the hexadecimal representation for a byte with all bits set to 1.
        self._sieve = bytearray([0xff] * sieve_size_bytes)

        # Helper function to clear a bit (mark a number as not prime).
        def _clear_bit(n):
            byte_index = n // 8
            bit_index = n % 8
            # Create a mask with the target bit as 0 and all others as 1, then AND it.
            mask = ~(1 << bit_index)
            self._sieve[byte_index] &= mask

        # Mark 0 and 1 as not prime.
        _clear_bit(0)
        _clear_bit(1)

        # Implement the Sieve of Eratosthenes.
        # We only need to iterate up to the square root of the maximum number.
        for p in range(2, int(math.sqrt(self._MAX_N)) + 1):
            # If p is prime (its bit is still 1), then mark its multiples as not prime.
            if self.isprime(p):
                # Start marking from p*p. Smaller multiples have already been marked
                # by smaller primes.
                for i in range(p * p, self._MAX_N, p):
                    _clear_bit(i)

    def isprime(self, p: int) -> bool:
        """
        Checks if p is a prime number in O(1) time.
        """
        if not (0 <= p < self._MAX_N):
            return False # Numbers outside the pre-computed range are not handled.
        
        byte_index = p // 8
        bit_index = p % 8
        mask = 1 << bit_index
        
        # Check if the bit is set using a bitwise AND.
        return (self._sieve[byte_index] & mask) != 0

    def primes(self, n: int) -> list[int]:
        """
        Returns a list of all prime numbers less than or equal to n in O(n) time.
        """
        if not (0 <= n < self._MAX_N):
            raise ValueError(f"Input n={n} is outside the pre-computed range of 0-{self._MAX_N-1}")
            
        # Iterate from 0 to n and collect all numbers for which isprime() is true.
        return [i for i in range(n + 1) if self.isprime(i)]

# --- Main execution block ---

# Create an instance of the data structure.
prime_checker = PrimeDataStructure(max_val=9999)

# 1. Demonstrate isprime(p)
p1 = 9973 # A large prime
p2 = 9974 # A non-prime
print(f"Is {p1} prime? {prime_checker.isprime(p1)}")
print(f"Is {p2} prime? {prime_checker.isprime(p2)}")
print("-" * 20)

# 2. Demonstrate primes(n)
n = 30
print(f"Primes up to {n}: {prime_checker.primes(n)}")
print("-" * 20)

# 3. Calculate and print the maximal size of the data structure.
# The maximal size is determined by the number of elements we need to check (0 to 9999).
NUM_ELEMENTS = 10000
BITS_PER_BYTE = 8
max_size_in_bytes = (NUM_ELEMENTS + BITS_PER_BYTE - 1) // BITS_PER_BYTE

print("Calculation for maximal size in bytes:")
print(f"The data structure needs to store primality for {NUM_ELEMENTS} numbers (0-9999).")
print(f"Using a bit array, each number requires 1 bit.")
print(f"The size in bytes is the ceiling of the total bits divided by {BITS_PER_BYTE}.")
print(f"Equation: ceil( {NUM_ELEMENTS} / {BITS_PER_BYTE} ) = {max_size_in_bytes}")
<<<1250>>>