import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations.
    It uses a pre-computed bit array (via Sieve of Eratosthenes) to store
    primality information for all numbers up to a given limit.
    """
    def __init__(self, limit=10000):
        """
        Initializes the data structure by building a prime sieve up to the limit.
        """
        if limit < 2:
            self.sieve = bytearray()
            self.limit = limit
            return

        self.limit = limit
        # A bit array to store primality. 1 = prime, 0 = composite.
        # Size in bytes is ceil(limit / 8).
        # We use (limit + 7) // 8 for integer division ceiling.
        self.sieve = bytearray([0b11111111] * ((limit + 7) // 8))

        # Helper function to set a bit to 0 (mark as not prime)
        def mark_not_prime(n):
            byte_index = n >> 3  # n // 8
            bit_index = n & 7    # n % 8
            self.sieve[byte_index] &= ~(1 << bit_index)

        # 0 and 1 are not prime.
        mark_not_prime(0)
        mark_not_prime(1)

        # Sieve of Eratosthenes
        # Iterate from 2 up to sqrt(limit)
        for p in range(2, int(math.sqrt(limit)) + 1):
            # If p is prime (its bit is still 1)
            if self.isprime(p):
                # Mark all multiples of p as not prime
                for multiple in range(p * p, limit, p):
                    mark_not_prime(multiple)

    def isprime(self, p: int) -> bool:
        """
        Checks if a number p is prime.
        Time Complexity: O(1)
        """
        if not 0 <= p < self.limit:
            return False

        byte_index = p >> 3  # p // 8
        bit_index = p & 7    # p % 8
        
        # Check if the bit at the given position is 1
        return (self.sieve[byte_index] >> bit_index) & 1 == 1

    def primes(self, n: int) -> list:
        """
        Returns a list of all prime numbers less than or equal to n.
        Time Complexity: O(n)
        """
        if n >= self.limit:
            n = self.limit - 1
            
        prime_list = []
        # Iterate from 2 up to n, checking primality for each number.
        for i in range(2, n + 1):
            if self.isprime(i):
                prime_list.append(i)
        return prime_list

# Main execution block
if __name__ == "__main__":
    MAX_NUM = 10000
    
    print("Initializing the data structure for primes up to 10000...")
    prime_checker = PrimeDataStructure(MAX_NUM)
    print("Initialization complete.")
    print("-" * 30)

    # --- Demonstrate isprime(p) ---
    print("Demonstrating isprime(p) [O(1) complexity]:")
    print(f"isprime(9973): {prime_checker.isprime(9973)}")
    print(f"isprime(9974): {prime_checker.isprime(9974)}")
    print("-" * 30)

    # --- Demonstrate primes(n) ---
    print("Demonstrating primes(n) [O(n) complexity]:")
    n_value = 40
    prime_list_up_to_n = prime_checker.primes(n_value)
    print(f"Primes up to {n_value}: {prime_list_up_to_n}")
    print("-" * 30)
    
    # --- Calculate and print the maximal size of the data structure ---
    print("Calculating the maximal size of the data structure in bytes:")
    # The size is determined by the bit array (sieve) needed to store
    # primality for numbers from 0 to MAX_NUM - 1. We need MAX_NUM bits.
    size_in_bytes = len(prime_checker.sieve)
    
    # The prompt requires printing the equation with the final numbers
    print(f"The calculation is based on storing {MAX_NUM} bits for primality.")
    print(f"Number of bytes = ceiling({MAX_NUM} / 8)")
    print(f"Equation: ({MAX_NUM} + 7) // 8 = {size_in_bytes}")
    print(f"The maximal size of this data structure is {size_in_bytes} bytes.")