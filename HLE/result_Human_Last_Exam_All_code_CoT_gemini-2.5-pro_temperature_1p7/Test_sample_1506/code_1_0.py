import math

class PrimeHandler:
    """
    A data structure for efficient prime number operations using a bit array sieve.
    This structure is optimized for memory by using one bit per number.
    """

    def __init__(self, max_limit):
        """
        Initializes the prime handler by pre-computing primes up to max_limit
        using the Sieve of Eratosthenes on a bit array.
        """
        if max_limit < 0:
            raise ValueError("Maximum limit must be non-negative.")
            
        self.max_limit = max_limit
        # Calculate the size of the bit array in bytes. Each byte stores 8 bits.
        # We add 7 to max_limit before integer division to find the ceiling.
        num_bytes = (max_limit + 7) // 8
        
        # The core data structure: a bytearray acting as a bit array.
        # Initialize all bits to 1 (representing 'is prime').
        self._sieve = bytearray([0b11111111] * num_bytes)

        # Helper function to clear a bit (set to 0), marking a number as not prime.
        def _clear_bit(n):
            """Sets the nth bit to 0."""
            if n < self.max_limit:
                byte_index = n // 8
                bit_index = n % 8
                # Create a mask like 11111110 to clear the bit at bit_index
                mask = ~(1 << bit_index)
                self._sieve[byte_index] &= mask

        # Mark 0 and 1 as not prime
        if max_limit > 0:
            _clear_bit(0)
        if max_limit > 1:
            _clear_bit(1)

        # Apply the Sieve of Eratosthenes
        # We only need to check for primes up to sqrt(max_limit)
        for p in range(2, int(math.sqrt(max_limit)) + 1):
            if self.isprime(p):  # If p is prime
                # Mark all multiples of p as not prime.
                # Start from p*p, as smaller multiples were already handled.
                for multiple in range(p * p, max_limit, p):
                    _clear_bit(multiple)

    def isprime(self, p):
        """
        Checks if p is a prime number.
        Time complexity: O(1)
        """
        if not 0 <= p < self.max_limit:
             raise ValueError(f"Input {p} is out of the supported range [0, {self.max_limit-1}]")

        byte_index = p // 8
        bit_index = p % 8
        # Check if the bit at the given position is 1
        return (self._sieve[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all primes p <= n.
        Time complexity: O(n)
        """
        if not 0 <= n < self.max_limit:
            raise ValueError(f"Input {n} is out of the supported range [0, {self.max_limit-1}]")
            
        # Iterate from 2 to n and collect primes using the pre-computed sieve
        return [i for i in range(2, n + 1) if self.isprime(i)]

    def get_data_structure_size_in_bytes(self):
        """
        Returns the size of the core data structure (the bit array).
        """
        return len(self._sieve)

# --- Main execution block ---

# Define the maximum number for which we need to store primality.
# Since p and n < 10000, the maximum index needed is 9999.
# This requires a structure of size 10000 to cover indices 0-9999.
MAX_NUM = 10000

# 1. Initialize the data structure. The Sieve is computed here.
prime_handler = PrimeHandler(MAX_NUM)

# 2. Calculate the maximal size of the data structure in bytes.
# The data structure is the sieve (bytearray). Its size is fixed.
max_size_in_bytes = prime_handler.get_data_structure_size_in_bytes()

# 3. Print the explanation and the final size calculation.
print(f"The most memory-efficient data structure is a bit array.")
print(f"To store primality for numbers from 0 to {MAX_NUM-1}, we need {MAX_NUM} bits.")
print("\nThe size calculation is as follows:")
print(f"Number of bytes = ceil({MAX_NUM} bits / 8 bits_per_byte)")
# Show the integer arithmetic to get the ceiling
c_part_1 = MAX_NUM // 8
c_part_2 = int(MAX_NUM % 8 > 0)
print(f"                = {MAX_NUM} // 8 + ({MAX_NUM} % 8 > 0)")
print(f"                = {c_part_1} + {c_part_2}")
print(f"                = {max_size_in_bytes}")

print(f"\nThe maximal size of this data structure is {max_size_in_bytes} bytes.")

# # --- Demonstration of operations (optional) ---
# print("\n--- Demonstrating isprime(p) ---")
# test_prime = 9973 # A prime number < 10000
# test_composite = 9975
# print(f"isprime({test_prime}): {prime_handler.isprime(test_prime)}")
# print(f"isprime({test_composite}): {prime_handler.isprime(test_composite)}")
# print("\n--- Demonstrating primes(n) ---")
# primes_list = prime_handler.primes(30)
# print(f"primes(30): {primes_list}")
<<<1250>>>