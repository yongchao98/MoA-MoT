import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations
    using a pre-computed bit array (Sieve of Eratosthenes).
    """

    def __init__(self, limit=10000):
        """
        Initializes the data structure by running a sieve up to the limit.
        """
        if limit < 0:
            raise ValueError("Limit cannot be negative.")
        self.limit = limit
        
        # Calculate size for bit array (1 bit per number)
        # We need (limit + 7) // 8 bytes.
        byte_size = (limit + 7) // 8
        # Initialize all bits to 1 (all numbers assumed prime)
        self._sieve = bytearray(b'\xff') * byte_size

        # Helper methods for bit manipulation
        def _clear_bit(p):
            # Sets the bit for number p to 0
            self._sieve[p >> 3] &= ~(1 << (p & 7))

        # 0 and 1 are not prime numbers
        if limit > 0:
            _clear_bit(0)
        if limit > 1:
            _clear_bit(1)

        # Sieve of Eratosthenes
        # Iterate from 2 up to sqrt(limit)
        for p in range(2, int(math.sqrt(limit)) + 1):
            # If p is prime (its bit is still 1)
            if self.isprime(p):
                # Mark all multiples of p as not prime
                # Start from p*p, as smaller multiples are already marked.
                for i in range(p * p, limit, p):
                    _clear_bit(i)

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be less than the initialized limit.
        """
        if not 0 <= p < self.limit:
            raise ValueError(f"Input {p} is out of the supported range [0, {self.limit-1}]")
            
        # Check the p-th bit in the sieve
        return (self._sieve[p >> 3] >> (p & 7)) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all prime numbers less than or equal to n in O(n) time.
        n must be less than the initialized limit.
        """
        if not 0 <= n < self.limit:
             raise ValueError(f"Input {n} is out of the supported range [0, {self.limit-1}]")

        prime_list = []
        for i in range(2, n + 1):
            if self.isprime(i):
                prime_list.append(i)
        return prime_list

def main():
    """
    Main function to demonstrate the PrimeDataStructure and calculate its size.
    """
    # Initialize the data structure for numbers up to 10000
    prime_checker = PrimeDataStructure(10000)

    # --- Demonstrate the operations ---
    print("--- Demonstrating Operations ---")
    p1 = 9973 # A prime number
    p2 = 9974 # Not a prime number
    print(f"Is {p1} prime? {prime_checker.isprime(p1)}")
    print(f"Is {p2} prime? {prime_checker.isprime(p2)}")

    n = 30
    print(f"Primes up to {n}: {prime_checker.primes(n)}")
    print("-" * 20)

    # --- Calculate and print the maximal size ---
    print("--- Maximal Size Calculation ---")
    num_numbers = 10000
    bits_per_byte = 8
    
    # In Python, len() of a bytearray returns its size in bytes.
    # This directly corresponds to our calculation.
    actual_size = len(prime_checker._sieve)
    
    print(f"The data structure needs to store primality for {num_numbers} numbers (0-9999).")
    print("Using a bit array, this requires 1 bit per number.")
    print(f"The total number of bytes is calculated by dividing the total bits by the number of bits in a byte.")
    print("Equation: Total Bytes = Number of numbers / Bits per byte")
    
    # Print each number in the final equation as requested
    print(f"{num_numbers} / {bits_per_byte} = {actual_size}")
    
    print(f"\nThe maximal size of this data structure's core storage is {actual_size} bytes.")


if __name__ == "__main__":
    main()

<<<1250>>>