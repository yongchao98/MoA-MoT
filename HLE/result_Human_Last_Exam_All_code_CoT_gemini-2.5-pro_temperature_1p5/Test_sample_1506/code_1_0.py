import math
import sys

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations using a bit array sieve.
    - isprime(p): O(1) check
    - primes(n): O(n) list generation
    """
    def __init__(self, max_n=10000):
        """Initializes the data structure by pre-computing primes up to max_n."""
        if max_n < 0:
            raise ValueError("max_n must be non-negative")
        # The structure needs to handle numbers < 10000, so the range is 0-9999.
        # This means we need an array of size 10000.
        self.limit = max_n
        
        # A bit array is the most memory-efficient storage for boolean flags.
        # We use a bytearray, where each byte stores 8 flags.
        num_bytes = (self.limit + 7) // 8
        self.bit_array = bytearray([0xFF] * num_bytes) # Initialize all bits to 1 (prime)

        # --- Use Sieve of Eratosthenes to populate the bit array ---
        # Helper to set a bit to 0 (not prime) or 1 (prime)
        def _set_bit(n, value):
            byte_index, bit_index = n >> 3, n & 7
            mask = 1 << bit_index
            if value:
                self.bit_array[byte_index] |= mask
            else:
                self.bit_array[byte_index] &= ~mask

        # 0 and 1 are not prime.
        if self.limit > 0: _set_bit(0, 0)
        if self.limit > 1: _set_bit(1, 0)

        # Iterate from 2 up to sqrt(limit)
        for p in range(2, int(math.sqrt(self.limit)) + 1):
            # If p is prime (its bit is still 1)
            if self.isprime(p):
                # Mark all multiples of p as not prime (set bit to 0).
                # Start from p*p, as smaller multiples have been marked already.
                for i in range(p * p, self.limit, p):
                    _set_bit(i, 0)

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be < 10000.
        """
        if not 0 <= p < self.limit:
            return False 
        
        byte_index = p >> 3 # Equivalent to p // 8
        bit_index = p & 7   # Equivalent to p % 8
        mask = 1 << bit_index
        return (self.bit_array[byte_index] & mask) != 0

    def primes(self, n):
        """
        Returns a list of all prime numbers p <= n in O(n) time.
        n must be < 10000.
        """
        if not 0 <= n < self.limit:
            # Cap n to the highest pre-computed value.
            n = self.limit -1

        prime_list = []
        # 2 is the first prime. Handle it separately to optimize the loop.
        if n >= 2:
             prime_list.append(2)
        # Check only odd numbers starting from 3 for other primes.
        for i in range(3, n + 1, 2):
            if self.isprime(i):
                prime_list.append(i)
        return prime_list
    
    def get_data_structure_size(self):
        """Returns the size of the core data structure in bytes."""
        return sys.getsizeof(self.bit_array)

if __name__ == "__main__":
    MAX_LIMIT = 10000  # For numbers p < 10000, we need indices 0 through 9999.
    
    # Create the data structure instance.
    prime_sieve = PrimeDataStructure(max_n=MAX_LIMIT)

    # Calculate the size of the data structure.
    # The data structure itself is a bytearray that holds all the primality info.
    num_bits = MAX_LIMIT
    num_bytes_data = (num_bits + 7) // 8
    total_size_bytes = prime_sieve.get_data_structure_size()
    overhead = total_size_bytes - num_bytes_data

    print("To store the primality of all numbers up to 9999, we need a bit array.")
    print(f"This requires one bit for each of the {num_bits} numbers.")
    print("\nFinal Equation for Memory Size:")
    print(f"Data size = ceil({num_bits} / 8) = {num_bytes_data} bytes")
    print(f"Python Object Overhead = {overhead} bytes (this may vary slightly by system architecture)")
    print(f"Maximal Size of Data Structure = {num_bytes_data} + {overhead} = {total_size_bytes} bytes")
    
    # Example usage (optional, uncomment to test)
    # print(f"\nIs 9999 a prime? {prime_sieve.isprime(9999)}")
    # print(f"Is 9973 a prime? {prime_sieve.isprime(9973)}")
    # print(f"Primes up to 30: {prime_sieve.primes(30)}")
    
    print(f"\n<<< {total_size_bytes} >>>")
