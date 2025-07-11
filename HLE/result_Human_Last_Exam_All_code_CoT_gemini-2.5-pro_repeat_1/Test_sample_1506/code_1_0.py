import math
import sys

class PrimeSieve:
    """
    A memory-efficient data structure for prime number operations using a bit array,
    which is the most efficient in terms of memory. The primality of numbers up to
    a given limit is pre-computed using the Sieve of Eratosthenes.
    """
    def __init__(self, max_limit=10000):
        """
        Initializes the sieve, pre-computing primes up to max_limit.
        - Time complexity of initialization: O(N log log N) where N is max_limit.
        - Space complexity: O(N) bits, which is O(N/8) bytes.
        """
        if max_limit < 0:
            raise ValueError("Maximum limit must be non-negative.")
        self._max_limit = max_limit
        
        # Calculate the size of the bit array in bytes. Each byte stores 8 flags.
        num_bytes = (max_limit + 7) // 8
        
        # The core data structure: a bytearray acting as a bit array.
        # It's initialized to all 1s, assuming all numbers are prime initially.
        self._bits = bytearray([0xff] * num_bytes)

        # Helper to clear a bit (i.e., mark a number as not prime)
        def _clear_bit(n):
            byte_index = n // 8
            bit_index = n % 8
            # The mask `~(1 << bit_index)` has a 0 at bit_index and 1s elsewhere.
            # ANDing with this mask clears the target bit, leaving others unchanged.
            self._bits[byte_index] &= ~(1 << bit_index)

        # Mark 0 and 1 as not prime.
        if max_limit > 0:
            _clear_bit(0)
            _clear_bit(1)
            
        # Sieve of Eratosthenes algorithm
        # Iterate from 2 up to the square root of the limit.
        for number in range(2, int(math.sqrt(max_limit)) + 1):
            # If the number is still marked as prime, then all its multiples are not.
            if self.isprime(number):
                # Start marking from number*number. Multiples smaller than that
                # (e.g., 2*number, 3*number) would have been marked by smaller primes.
                for multiple in range(number * number, max_limit, number):
                    _clear_bit(multiple)

    def isprime(self, p):
        """
        Checks if a number p is prime.
        Time complexity: O(1).
        """
        if not 0 <= p < self._max_limit:
            raise ValueError(f"Input {p} is out of the pre-computed range [0, {self._max_limit-1}]")
            
        byte_index = p // 8
        bit_index = p % 8
        # Check if the corresponding bit is set to 1.
        return (self._bits[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all prime numbers less than or equal to n.
        Time complexity: O(n).
        """
        if not 0 <= n < self._max_limit:
            raise ValueError(f"Input {n} is out of the pre-computed range [0, {self._max_limit-1}]")
        
        # Iterate from 0 to n, collecting numbers that are prime using the O(1) isprime check.
        return [i for i in range(n + 1) if self.isprime(i)]
        
    def get_data_structure_size_in_bytes(self):
        """
        Returns the size of the core data structure (the bit array) in bytes.
        """
        return len(self._bits)

# --- Main execution block ---
# Create an instance of the data structure for numbers < 10000.
MAX_LIMIT = 10000
sieve = PrimeSieve(MAX_LIMIT)

# --- Example Usage ---
print("--- Example Operations ---")
p1 = 9973 # The largest prime < 10000
p2 = 9974
print(f"isprime({p1}): {sieve.isprime(p1)}")
print(f"isprime({p2}): {sieve.isprime(p2)}")

n = 30
print(f"primes({n}): {sieve.primes(n)}")


# --- Maximal Size Calculation ---
print("\n--- Data Structure Size Calculation ---")
size_in_bytes = sieve.get_data_structure_size_in_bytes()
print(f"The data structure needs to store primality for {MAX_LIMIT} numbers (from 0 to {MAX_LIMIT - 1}).")
print("The most memory-efficient method is a bit array, using 1 bit per number.")
print(f"Total bits required: {MAX_LIMIT}")
print(f"Number of bits in a byte: 8")
print("To find the number of bytes, we calculate: ceil(Total bits / Bits per byte)")
print(f"Calculation: ceil({MAX_LIMIT} / 8) = ({MAX_LIMIT} + 7) // 8 = {size_in_bytes}")
print(f"The maximal size of this data structure is {size_in_bytes} bytes.")

<<<1250>>>