import math
import sys

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations
    using a pre-computed sieve.
    """
    def __init__(self, max_num=10000):
        # The sieve needs to hold values up to max_num - 1.
        # We create it of size max_num to have indices 0 to max_num-1.
        self.max_num = max_num
        # Using a bytearray for memory efficiency compared to a standard list.
        # 1 represents a potential prime, 0 represents not a prime.
        self.sieve = bytearray([1]) * max_num
        # 0 and 1 are not prime numbers.
        if max_num > 1:
            self.sieve[0] = 0
            self.sieve[1] = 0

        # Sieve of Eratosthenes algorithm to populate the sieve.
        # We only need to iterate up to the square root of max_num.
        for i in range(2, int(math.sqrt(max_num))):
            if self.sieve[i] == 1:
                # Mark all multiples of i as not prime.
                # Start marking from i*i.
                for j in range(i * i, max_num, i):
                    self.sieve[j] = 0

    def isprime(self, p):
        """
        Checks if a number p is prime.
        Time Complexity: O(1)
        """
        if 0 <= p < self.max_num:
            return self.sieve[p] == 1
        return False # Out of pre-computed bounds

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n.
        Time Complexity: O(n)
        """
        if n >= self.max_num:
            n = self.max_num -1

        prime_list = []
        for i in range(n + 1):
            if self.sieve[i] == 1:
                prime_list.append(i)
        return prime_list

def main():
    # Instantiate the data structure
    prime_handler = PrimeDataStructure(max_num=10000)

    # --- Demonstrate the operations ---
    # 1. isprime(p)
    p = 9973 # A large prime
    print(f"Is {p} a prime? {prime_handler.isprime(p)}")
    p = 9975 # Not a prime
    print(f"Is {p} a prime? {prime_handler.isprime(p)}\n")

    # 2. primes(n)
    n = 100
    print(f"Primes up to {n}: {prime_handler.primes(n)}\n")


    # --- Calculate the maximal size of the most efficient data structure ---
    print("Calculating the size of the most efficient data structure (a bit array):")
    # The structure needs to store a boolean flag for each number from 0 to 9999.
    flags_needed = 10000
    # The most efficient way to store a boolean flag is a single bit.
    # There are 8 bits in a byte.
    bits_per_byte = 8
    # The size in bytes is the number of flags divided by bits per byte.
    # We use ceiling division in case the number of flags is not a multiple of 8.
    size_in_bytes = math.ceil(flags_needed / bits_per_byte)

    print(f"The data structure needs to store {flags_needed} boolean flags (one for each number from 0 to 9999).")
    print(f"The most memory-efficient representation uses 1 bit per flag.")
    print(f"There are {bits_per_byte} bits in a byte.")
    print(f"Calculation for size in bytes: {flags_needed} / {bits_per_byte} = {size_in_bytes}")
    print(f"Maximal size of this data structure: {size_in_bytes} bytes.")

if __name__ == "__main__":
    main()

<<<1250>>>