import math

class PrimeChecker:
    """
    Designs and implements a memory-efficient data structure for prime number operations.
    The core data structure is a boolean list representing a bit array,
    pre-populated using the Sieve of Eratosthenes.
    """

    def __init__(self, max_num=10000):
        """
        Initializes the primality lookup table for numbers up to max_num.
        """
        if max_num < 0:
            raise ValueError("max_num must be non-negative")
            
        self._max_num = max_num
        # Create a boolean list. True means the index is potentially prime.
        self._is_prime = [True] * self._max_num
        
        # 0 and 1 are not prime numbers
        if self._max_num > 0:
            self._is_prime[0] = False
        if self._max_num > 1:
            self._is_prime[1] = False

        # Sieve of Eratosthenes algorithm to mark non-primes
        for i in range(2, int(math.sqrt(self._max_num)) + 1):
            if self._is_prime[i]:
                # Mark all multiples of i as not prime
                for multiple in range(i * i, self._max_num, i):
                    self._is_prime[multiple] = False

    def isprime(self, p: int) -> bool:
        """
        Checks if p is a prime number in O(1) time.
        Args:
            p: An integer where 0 <= p < 10000.
        Returns:
            True if p is prime, False otherwise.
        """
        if not 0 <= p < self._max_num:
            raise ValueError(f"Input p must be between 0 and {self._max_num - 1}")
        return self._is_prime[p]

    def primes(self, n: int) -> list[int]:
        """
        Returns a list of all primes p <= n in O(n) time.
        Args:
            n: An integer where 0 <= n < 10000.
        Returns:
            A list of prime numbers.
        """
        if not 0 <= n < self._max_num:
            raise ValueError(f"Input n must be between 0 and {self._max_num - 1}")
        
        # Iterate from 2 to n and collect primes using the pre-computed table
        return [i for i in range(2, n + 1) if self._is_prime[i]]
        
    def get_data_structure_size_in_bytes(self) -> int:
        """
        Calculates the size of the most memory-efficient representation (bit array)
        of the primality lookup table in bytes.
        """
        # The lookup table requires one bit for each number from 0 to _max_num - 1.
        num_bits = self._max_num
        # To get the size in bytes, divide by 8 and take the ceiling.
        num_bytes = math.ceil(num_bits / 8)
        return int(num_bytes)

def main():
    """
    Main function to demonstrate PrimeChecker and calculate its size.
    """
    # Initialize the data structure for numbers < 10000
    limit = 10000
    prime_checker = PrimeChecker(max_num=limit)
    
    # --- Demonstrate the functionality ---
    
    # 1. isprime(p)
    p1 = 9973  # A prime number
    p2 = 9974  # Not a prime number
    print(f"Is the number {p1} a prime? {prime_checker.isprime(p1)}")
    print(f"Is the number {p2} a prime? {prime_checker.isprime(p2)}")
    print("-" * 20)

    # 2. primes(n)
    n = 30
    prime_list = prime_checker.primes(n)
    print(f"The list of prime numbers up to {n} is: {prime_list}")
    print("-" * 20)
    
    # --- Calculate and print the maximal size ---
    
    print("Calculating the maximal size of the data structure:")
    size_in_bytes = prime_checker.get_data_structure_size_in_bytes()
    num_bits = prime_checker._max_num
    bits_per_byte = 8
    
    print(f"The most efficient data structure is a bit array.")
    print(f"It needs to store a boolean flag for each number from 0 to {num_bits - 1}.")
    print(f"Total bits required: {num_bits}")
    print(f"There are {bits_per_byte} bits in a byte.")
    print(f"Equation for size in bytes: ceil({num_bits} / {bits_per_byte})")
    print(f"Calculated size: {size_in_bytes} bytes")


if __name__ == "__main__":
    main()
