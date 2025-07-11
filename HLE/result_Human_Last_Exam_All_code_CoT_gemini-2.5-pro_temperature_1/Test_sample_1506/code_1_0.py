import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations using a
    pre-computed bit array filled by the Sieve of Eratosthenes.
    """

    def __init__(self, limit=10000):
        """
        Initializes the data structure by running the Sieve of Eratosthenes
        to populate a bit array.
        """
        if limit < 0:
            raise ValueError("Limit must be a non-negative number.")
        self.limit = limit
        # Each byte stores 8 numbers' primality status.
        # We use ceiling division: (limit + 7) // 8
        num_bytes = (limit + 7) // 8
        
        # The core data structure: a bit array.
        # Initialize all bits to 1 (True), assuming all are prime initially.
        self._bit_array = bytearray([0xFF] * num_bytes)

        # Helper method to clear a bit (mark a number as not prime)
        def _clear_bit(n):
            byte_index = n // 8
            bit_index = n % 8
            # The '&=' operation with a bitwise NOT '~' clears the specific bit.
            self._bit_array[byte_index] &= ~(1 << bit_index)

        # 0 and 1 are not prime numbers.
        if limit > 0:
            _clear_bit(0)
        if limit > 1:
            _clear_bit(1)

        # Sieve of Eratosthenes algorithm
        # We only need to sieve up to the square root of the limit.
        for p in range(2, int(math.sqrt(limit)) + 1):
            if self.isprime(p):  # If p is still marked as prime
                # Mark all multiples of p (starting from p*p) as not prime.
                for multiple in range(p * p, limit, p):
                    _clear_bit(multiple)

    def isprime(self, p):
        """
        Checks if p is a prime number in O(1) time.
        p must be less than the initialized limit.
        """
        if not 0 <= p < self.limit:
            # For this problem, we assume valid inputs as per the prompt.
            # In a real-world scenario, robust error handling is needed.
            return False
            
        # Check the p-th bit in the bit array.
        byte_index = p // 8
        bit_index = p % 8
        return (self._bit_array[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n in O(n) time.
        n must be less than the initialized limit.
        """
        if not 0 <= n < self.limit:
             return []
        
        # isprime() is O(1), so a list comprehension iterating n times is O(n).
        return [i for i in range(2, n + 1) if self.isprime(i)]


# --- Main execution block to demonstrate usage and calculate size ---

# Per the problem, the maximum value for p and n is less than 10000.
MAX_LIMIT = 10000

# 1. Initialize the data structure.
prime_db = PrimeDataStructure(limit=MAX_LIMIT)

# 2. Demonstrate the isprime(p) operation.
print(f"Is 9973 a prime number? {prime_db.isprime(9973)}")
print(f"Is 100 a prime number? {prime_db.isprime(100)}")
print("-" * 20)

# 3. Demonstrate the primes(n) operation.
print(f"Primes up to 30: {prime_db.primes(30)}")
print("-" * 20)

# 4. Calculate and print the maximal size of the data structure.
print("Calculating the maximal size of the data structure:")
print("The structure must store the primality for numbers from 0 to 9999.")
print("This is a total of 10000 numbers.")
print("Using a bit array, we store one bit for each number.")
print("Total bits required = 10000")
print("To convert bits to bytes, we divide by 8 (since 1 byte = 8 bits).")

bits_needed = 10000
bits_per_byte = 8
bytes_needed = bits_needed // bits_per_byte

print(f"The calculation is: {bits_needed} / {bits_per_byte} = {bytes_needed}")
print(f"So, the maximal size of this data structure is {bytes_needed} bytes.")

<<<1250>>>