import sys

class PrimeChecker:
    """
    A memory-efficient data structure for prime number operations using a pre-computed bit sieve.

    This class pre-computes the primality of all numbers up to a specified maximum
    limit using a Sieve of Eratosthenes stored in a memory-efficient bit array
    (implemented with a bytearray).

    Operations:
    - isprime(p): Checks if a number p is prime in O(1) time.
    - primes(n): Returns a list of all primes up to n in O(n) time.
    """
    def __init__(self, max_num=10000):
        """
        Initializes the sieve for all numbers from 0 to max_num-1.
        """
        if max_num <= 1:
            raise ValueError("max_num must be greater than 1")
        self.max_limit = max_num
        # Calculate the size of the bit array in bytes. (max_num + 7) // 8 rounds up.
        self.sieve_size_bytes = (self.max_limit + 7) // 8
        # Initialize bytearray with all bits set to 1 (representing potentially prime)
        self.sieve = bytearray([0xff] * self.sieve_size_bytes)

        # Mark 0 and 1 as not prime by clearing their corresponding bits
        self._clear_bit(0)
        self._clear_bit(1)

        # Pre-compute primes using Sieve of Eratosthenes
        # We only need to sieve up to the square root of the limit
        for i in range(2, int(self.max_limit**0.5) + 1):
            if self.isprime(i):
                # Mark all multiples of i (starting from i*i) as not prime
                for multiple in range(i * i, self.max_limit, i):
                    self._clear_bit(multiple)

    def _clear_bit(self, p):
        """Helper function to clear the bit for number p (marks it as not prime)."""
        byte_index = p // 8
        bit_index = p % 8
        # Create a mask to clear only the specific bit (e.g., ~(0b00010000))
        # and apply it with a bitwise AND
        mask = ~(1 << bit_index)
        self.sieve[byte_index] &= mask

    def isprime(self, p):
        """
        Checks if p is a prime number in O(1) time.
        p must be less than the max_limit used for initialization (10000).
        """
        if not 0 <= p < self.max_limit:
            raise ValueError(f"Input {p} is out of the pre-computed range [0, {self.max_limit-1}]")

        byte_index = p // 8
        bit_index = p % 8
        # Check if the specific bit is set to 1
        return (self.sieve[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all prime numbers p <= n in O(n) time.
        n must be less than the max_limit used for initialization (10000).
        """
        if not 0 <= n < self.max_limit:
            raise ValueError(f"Input {n} is out of the pre-computed range [0, {self.max_limit-1}]")
        
        # This list comprehension iterates from 0 to n, performing an O(1) check
        # for each number, resulting in O(n) total time complexity.
        return [i for i in range(n + 1) if self.isprime(i)]

    def get_data_structure_size(self):
        """
        Calculates and returns the size of the main data structure (the sieve) in bytes.
        """
        return sys.getsizeof(self.sieve)


def main():
    """
    Main function to instantiate the PrimeChecker, demonstrate its use,
    and calculate the maximal size of its data structure.
    """
    # Initialize the data structure for primes up to 9999 (max_num=10000)
    max_number = 10000
    prime_checker = PrimeChecker(max_num=max_number)

    print("--- Demonstrating Operations ---")
    print(f"Is 9973 prime? {prime_checker.isprime(9973)}")
    print(f"Is 100 prime? {prime_checker.isprime(100)}")
    print(f"Primes up to 30: {prime_checker.primes(30)}")
    print("-" * 20)

    # --- Calculate and print the maximal size ---
    print("\n--- Calculating Maximal Size of the Data Structure ---")
    
    # The data structure is a bytearray that stores a bit for each number.
    # The equation for the number of bytes needed for the data is: ceil(number_of_bits / 8)
    # In integer arithmetic, this is calculated as (number_of_bits + 7) // 8.
    num_bits = max_number
    num_bytes_for_data = (num_bits + 7) // 8

    print(f"1. Total numbers to represent: {num_bits} (from 0 to 9999)")
    print(f"2. Bits per number: 1")
    print(f"3. Bits per byte: 8")
    print(f"4. Equation for minimum data bytes: ceil({num_bits} / {8}) = ({num_bits} + {8-1}) // {8} = {num_bytes_for_data} bytes")

    # sys.getsizeof() includes the Python object's overhead in addition to the data size.
    actual_size_in_bytes = prime_checker.get_data_structure_size()
    overhead = actual_size_in_bytes - num_bytes_for_data

    print(f"5. The bytearray data requires {num_bytes_for_data} bytes.")
    print(f"6. The Python bytearray object has an additional overhead of {overhead} bytes.")
    print("\n   Maximal size = (Data Size) + (Object Overhead)")
    print(f"   Maximal size = {num_bytes_for_data} + {overhead} = {actual_size_in_bytes} bytes")
    
    print(f"<<<{actual_size_in_bytes}>>>")


if __name__ == "__main__":
    main()