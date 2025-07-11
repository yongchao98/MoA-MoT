import math

class PrimeChecker:
    """
    A memory-efficient data structure for prime number operations using a bit array.
    """
    def __init__(self, max_val=9999):
        """
        Initializes the data structure, pre-computing primes up to max_val
        using the Sieve of Eratosthenes on a bit array.
        """
        self.max_val = max_val
        # Each byte holds 8 bits. We need to store max_val+1 numbers (0 to max_val).
        # The size of the bytearray is ceil((max_val + 1) / 8).
        size_in_bytes = (self.max_val + 8) // 8
        
        # This bit array is our core data structure.
        # Initialize all bits to 1 (assume all numbers are prime initially).
        self.sieve = bytearray([0b11111111] * size_in_bytes)

        # Helper function to set a bit to 0 (mark as not prime)
        def _clear_bit(n):
            byte_index = n // 8
            bit_pos = n % 8
            # Create a mask to clear the specific bit, e.g., 11110111
            mask = ~(1 << bit_pos)
            self.sieve[byte_index] &= mask

        # 0 and 1 are not prime.
        _clear_bit(0)
        _clear_bit(1)

        # Sieve of Eratosthenes
        for i in range(2, int(math.sqrt(self.max_val)) + 1):
            if self.isprime(i):
                # Mark all multiples of i (starting from i*i) as not prime.
                for j in range(i * i, self.max_val + 1, i):
                    _clear_bit(j)

    def isprime(self, p: int) -> bool:
        """
        Checks if a number p is prime in O(1) time.
        p must be less than 10000.
        """
        if not (0 <= p <= self.max_val):
            raise ValueError(f"Input p must be between 0 and {self.max_val}")
        
        byte_index = p // 8
        bit_pos = p % 8
        # Check if the bit at the given position is 1.
        return (self.sieve[byte_index] >> bit_pos) & 1 == 1

    def primes(self, n: int) -> list[int]:
        """
        Returns a list of all prime numbers p <= n in O(n) time.
        n must be less than 10000.
        """
        if not (0 <= n <= self.max_val):
            raise ValueError(f"Input n must be between 0 and {self.max_val}")
        
        # Iterate from 2 to n and collect numbers if they are prime.
        return [i for i in range(2, n + 1) if self.isprime(i)]

    def get_data_structure_size(self) -> int:
        """
        Returns the size of the core data structure (the sieve) in bytes.
        """
        return len(self.sieve)

if __name__ == '__main__':
    # Initialize the prime checker. The problem states p, n < 10000,
    # so the structure must handle numbers up to 9999.
    prime_checker = PrimeChecker(max_val=9999)

    # --- Demonstrate isprime(p) ---
    print(f"isprime(9998): {prime_checker.isprime(9998)}") # Should be False
    print(f"isprime(9973): {prime_checker.isprime(9973)}") # Should be True
    print("-" * 20)

    # --- Demonstrate primes(n) ---
    list_of_primes = prime_checker.primes(n=100)
    print(f"Primes up to 100: {list_of_primes}")
    print("-" * 20)

    # --- Calculate and print the maximal size of the data structure ---
    print("The data structure needs to store primality for 10000 numbers (0 to 9999).")
    print("Using a bit array, each number requires 1 bit.")
    print("There are 8 bits in a byte.")
    
    num_of_numbers = 10000
    bits_in_byte = 8
    total_size_in_bytes = prime_checker.get_data_structure_size()
    
    # Print the equation as requested
    print(f"Maximal size calculation: {num_of_numbers} / {bits_in_byte} = {total_size_in_bytes}")

<<<1250>>>