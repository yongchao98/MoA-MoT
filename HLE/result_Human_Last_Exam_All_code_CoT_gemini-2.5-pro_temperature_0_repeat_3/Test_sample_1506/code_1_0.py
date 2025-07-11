import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations using a bit array.
    - isprime(p): O(1)
    - primes(n): O(n)
    """
    def __init__(self, max_num=9999):
        """
        Initializes the data structure by pre-computing primes up to max_num
        using the Sieve of Eratosthenes on a bit array.
        """
        self.max_num = max_num
        # Each byte holds 8 bits (flags for 8 numbers).
        # The size is ceiling of (max_num + 1) / 8.
        self.size_in_bytes = (max_num + 1 + 7) // 8
        # Initialize bytearray with all bits set to 1 (all numbers assumed prime).
        # 0xff is 11111111 in binary.
        self.bit_array = bytearray([0xff] * self.size_in_bytes)

        # Helper functions to set/get bits
        def clear_bit(n):
            # Sets the bit for number n to 0 (not prime)
            byte_index = n // 8
            bit_index = n % 8
            self.bit_array[byte_index] &= ~(1 << bit_index)

        def is_bit_set(n):
            # Checks if the bit for number n is 1 (is prime)
            byte_index = n // 8
            bit_index = n % 8
            return (self.bit_array[byte_index] >> bit_index) & 1

        # 0 and 1 are not prime
        clear_bit(0)
        clear_bit(1)

        # Sieve of Eratosthenes
        for number in range(2, int(math.sqrt(max_num)) + 1):
            if is_bit_set(number):  # If number is prime
                # Mark all multiples of the number as not prime
                for multiple in range(number * number, max_num + 1, number):
                    clear_bit(multiple)

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be within the pre-computed range [0, max_num].
        """
        if not 0 <= p <= self.max_num:
            raise ValueError(f"Input {p} is out of the supported range [0, {self.max_num}]")
        
        byte_index = p // 8
        bit_index = p % 8
        return (self.bit_array[byte_index] >> bit_index) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all prime numbers less than or equal to n in O(n) time.
        n must be within the pre-computed range [0, max_num].
        """
        if not 0 <= n <= self.max_num:
            raise ValueError(f"Input {n} is out of the supported range [0, {self.max_num}]")
        
        prime_list = []
        for i in range(n + 1):
            if self.isprime(i):
                prime_list.append(i)
        return prime_list

# --- Main execution ---
if __name__ == "__main__":
    # Initialize the data structure
    prime_checker = PrimeDataStructure(max_num=9999)

    # 1. Demonstrate isprime(p)
    p = 9973 # A prime number
    print(f"Is {p} a prime? {prime_checker.isprime(p)}")
    p = 9974 # Not a prime number
    print(f"Is {p} a prime? {prime_checker.isprime(p)}")
    print("-" * 20)

    # 2. Demonstrate primes(n)
    n = 30
    prime_list_up_to_n = prime_checker.primes(n)
    print(f"Primes up to {n}: {prime_list_up_to_n}")
    print("-" * 20)

    # 3. Calculate and print the maximal size of the data structure
    max_size = prime_checker.size_in_bytes
    print(f"The data structure needs to store primality for 10000 numbers (0 to 9999).")
    print(f"Using a bit array, each number requires 1 bit.")
    print(f"Total bits required = 10000")
    print(f"Total bytes required = 10000 / 8")
    print(f"Maximal size of the data structure is {max_size} bytes.")

    # Final answer in the required format
    print("<<<1250>>>")