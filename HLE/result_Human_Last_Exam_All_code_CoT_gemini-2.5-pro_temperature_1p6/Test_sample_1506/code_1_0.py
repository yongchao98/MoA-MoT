import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations up to a limit.
    
    It uses a bit array (implemented with a bytearray) which is pre-computed
    using the Sieve of Eratosthenes. This allows for O(1) primality tests
    and O(n) generation of primes up to n.
    """

    def __init__(self, limit=10000):
        """
        Initializes the data structure, running the Sieve of Eratosthenes.
        """
        if limit < 0:
            raise ValueError("Limit must be a non-negative integer.")
        
        # The limit of the numbers we can check (0 to limit-1)
        self.limit = limit
        
        # Calculate the size of the bytearray needed to store the bits.
        # Each byte stores 8 numbers' primality.
        self.size_in_bytes = (self.limit + 7) // 8
        
        # The bit array. Initialize all bits to 1 (True/prime).
        # 0xFF is a byte where all 8 bits are 1s.
        self._sieve = bytearray([0xFF] * self.size_in_bytes)

        # Mark 0 and 1 as not prime by clearing their corresponding bits.
        self._clear_bit(0)
        self._clear_bit(1)
        
        # Use Sieve of Eratosthenes to mark non-prime numbers
        for p in range(2, int(math.sqrt(limit)) + 1):
            if self._is_bit_set(p): # If p's bit is still 1, it's a prime
                # Mark all multiples of p as not prime (set their bits to 0)
                # Start marking from p*p
                for i in range(p * p, self.limit, p):
                    self._clear_bit(i)

    def _is_bit_set(self, n):
        """Checks if the bit corresponding to number n is 1."""
        # Find which byte the bit is in (n // 8)
        byte_index = n >> 3
        # Find the position of the bit within that byte (n % 8)
        bit_index = n & 7
        # Return True if the bit is 1, False otherwise
        return (self._sieve[byte_index] >> bit_index) & 1

    def _clear_bit(self, n):
        """Sets the bit corresponding to number n to 0."""
        byte_index = n >> 3
        bit_index = n & 7
        # Create a mask to clear only the specific bit
        # e.g., for bit 2, mask is ~(1 << 2) = ~4 = ~00000100 = 11111011
        mask = ~(1 << bit_index)
        self._sieve[byte_index] &= mask

    def isprime(self, p):
        """
        Checks if p is a prime number in O(1) time.
        
        Args:
            p: An integer where 0 <= p < self.limit.
        
        Returns:
            True if p is prime, False otherwise.
        """
        if not 0 <= p < self.limit:
            raise ValueError(f"Input p must be between 0 and {self.limit - 1}")
        return self._is_bit_set(p) == 1

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n in O(n) time.
        
        Args:
            n: An integer where 0 <= n < self.limit.
            
        Returns:
            A list of prime numbers.
        """
        if not 0 <= n < self.limit:
            raise ValueError(f"Input n must be between 0 and {self.limit - 1}")
        
        # Iterate from 2 to n and collect all numbers that are prime
        return [i for i in range(n + 1) if self.isprime(i)]

    def get_data_structure_size(self):
        """Returns the size of the core _sieve data structure in bytes."""
        return self.size_in_bytes

if __name__ == '__main__':
    # The problem specifies a limit for p and n up to 10000.
    MAX_LIMIT = 10000
    
    # This will pre-calculate all primes up to 9999
    prime_checker = PrimeDataStructure(MAX_LIMIT)

    # Calculate and print the maximal size of the data structure.
    # The question requires printing the numbers in the final equation.
    size = prime_checker.get_data_structure_size()
    print("Calculation for the maximal size of the data structure in bytes:")
    print(f"Numbers to check: 0 to {MAX_LIMIT - 1} ({MAX_LIMIT} numbers)")
    print(f"Bits needed: {MAX_LIMIT}")
    print(f"Bytes needed = ceil(Bits / 8) = ceil({MAX_LIMIT} / 8)")
    print(f"The equation in integer arithmetic is ({MAX_LIMIT} + 7) // 8")
    print(f"Size = ({MAX_LIMIT} + 7) // 8 = {size} bytes.")
    
    # Demonstrate that the methods work correctly.
    print("\n--- Functionality Demonstration ---")
    p_test = 9973 # A large prime
    n_test = 40
    print(f"isprime({p_test}): {prime_checker.isprime(p_test)}")
    print(f"isprime({p_test+1}): {prime_checker.isprime(p_test+1)}")
    print(f"primes({n_test}): {prime_checker.primes(n_test)}")
