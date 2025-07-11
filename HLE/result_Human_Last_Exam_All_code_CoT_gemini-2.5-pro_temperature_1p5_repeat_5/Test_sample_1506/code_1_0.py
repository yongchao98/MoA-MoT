import math

class PrimeSieve:
    """
    A data structure for efficient prime number operations using a pre-computed sieve.
    This implementation uses a bytearray for simplicity, where 1 byte is used per number.
    The most memory-efficient design, a bit array, would use only 1 bit per number.
    """
    def __init__(self, limit=10000):
        """
        Initializes the sieve for all numbers up to the limit using the
        Sieve of Eratosthenes algorithm.
        """
        if limit < 2:
            raise ValueError("Limit must be at least 2.")
        self.limit = limit
        # 1 represents a prime, 0 represents a non-prime.
        self._sieve = bytearray([1]) * limit
        self._sieve[0] = 0
        self._sieve[1] = 0
        for number in range(2, int(math.sqrt(limit))):
            if self._sieve[number] == 1:
                # Mark all multiples of this prime as not prime.
                # Start marking from number*number.
                for multiple in range(number * number, limit, number):
                    self._sieve[multiple] = 0

    def isprime(self, p):
        """
        Checks if a number is prime in O(1) time.
        
        Args:
            p: The number to check.
        
        Returns:
            True if p is prime, False otherwise.
        """
        if not 0 <= p < self.limit:
            raise ValueError(f"Input {p} is out of the supported range [0, {self.limit-1}]")
        return self._sieve[p] == 1

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n in O(n) time.
        
        Args:
            n: The upper bound (inclusive).
        
        Returns:
            A list of prime numbers up to n.
        """
        if not 0 <= n < self.limit:
             raise ValueError(f"Input {n} is out of the supported range [0, {self.limit-1}]")
        return [i for i in range(2, n + 1) if self._sieve[i] == 1]


# Step 1: Demonstrate the data structure's functionality.
prime_checker = PrimeSieve(10000)

print("Demonstrating the data structure operations:")
print(f"isprime(9973): {prime_checker.isprime(9973)}")
print(f"isprime(100): {prime_checker.isprime(100)}")
print(f"primes(30): {prime_checker.primes(30)}")


# Step 2: Calculate the maximal size of the most efficient data structure (bit array).
print("\nCalculating the size of the most memory-efficient data structure:")
max_numbers = 10000  # For numbers 0 to 9999
bits_per_number = 1
total_bits = max_numbers * bits_per_number
bits_in_a_byte = 8
total_bytes = total_bits / bits_in_a_byte

print(f"To store primality for {max_numbers} numbers, we need one bit per number.")
print(f"Total bits = {max_numbers} numbers * {bits_per_number} bit/number = {total_bits} bits")
print(f"Total bytes = {total_bits} bits / {bits_in_a_byte} bits/byte = {int(total_bytes)} bytes")
<<<1250>>>