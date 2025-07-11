import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for primality tests using a pre-computed bit sieve.
    It satisfies:
    - isprime(p): O(1)
    - primes(n): O(n)
    """
    def __init__(self, max_num=10000):
        """
        Initializes the data structure, pre-computing primes up to max_num
        using a Sieve of Eratosthenes on a bit array.
        """
        if max_num < 0:
            raise ValueError("max_num must be non-negative")
            
        self.max_num = max_num
        # Size of the bit array in bytes. (max_num + 7) // 8 performs ceiling division.
        self.size_in_bytes = (self.max_num + 7) // 8
        # Use a bytearray as a bit array, initialized with all bits as 1 (prime).
        # 0xFF is hexadecimal for 0b11111111.
        self.sieve = bytearray([0xFF] * self.size_in_bytes)

        # Helper to set a bit for a given number to 0 (not prime).
        def _clear_bit(n):
            byte_index = n // 8
            bit_index = n % 8
            # The & operator with the inverted mask clears the bit.
            self.sieve[byte_index] &= ~(1 << bit_index)

        # 0 and 1 are not prime.
        if self.max_num > 0:
            _clear_bit(0)
        if self.max_num > 1:
            _clear_bit(1)
        
        # Sieve of Eratosthenes algorithm
        # We only need to sieve up to the square root of max_num.
        for p in range(2, int(math.sqrt(self.max_num))):
            if self.isprime(p):
                # Mark all multiples of p as not prime.
                # Start marking from p*p, as smaller multiples are already marked.
                for i in range(p * p, self.max_num, p):
                    _clear_bit(i)

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be within the pre-computed range [0, self.max_num - 1].
        """
        if not 0 <= p < self.max_num:
            return False # Numbers outside the pre-computed range are not handled.
            
        byte_index = p // 8
        bit_index = p % 8
        # Check if the bit at the position for p is 1.
        return (self.sieve[byte_index] & (1 << bit_index)) != 0

    def primes(self, n):
        """
        Returns a list of all primes p such that p <= n in O(n) time.
        n must be within the pre-computed range [0, self.max_num - 1].
        """
        if n >= self.max_num:
            # Handle requests for n larger than our capacity by capping it.
            n = self.max_num - 1
            
        # List comprehension iterates from 2 to n, checking primality for each.
        return [p for p in range(2, n + 1) if self.isprime(p)]

# Main part of the script
if __name__ == '__main__':
    # The upper limit for primes is less than 10000.
    MAX_LIMIT = 10000

    # You can uncomment the following lines to test the functionality.
    # print("Initializing the data structure...")
    # prime_checker = PrimeDataStructure(MAX_LIMIT)
    # print("Initialization complete.")
    # print(f"\nIs 9973 prime? {prime_checker.isprime(9973)}")
    # print(f"Is 9975 prime? {prime_checker.isprime(9975)}")
    # print(f"Primes up to 30: {prime_checker.primes(30)}\n")

    # --- Calculate and print the maximal size of the data structure ---
    print("Calculating the maximal size of the data structure in bytes:")
    
    # The data structure needs to handle numbers p < 10000, so from 0 to 9999.
    total_numbers = 10000
    
    # Using a bit array, each flag takes 1 bit.
    bits_per_number = 1
    total_bits_required = total_numbers * bits_per_number
    
    # There are 8 bits in a byte.
    bits_in_a_byte = 8
    
    # To find the number of bytes, we divide the total bits by 8 and take the ceiling.
    size_in_bytes = (total_bits_required + bits_in_a_byte - 1) // bits_in_a_byte

    print(f"\nTotal numbers to store (0 to {total_numbers - 1}): {total_numbers}")
    print(f"Bits required per number: {bits_per_number}")
    print(f"Bits in one byte: {bits_in_a_byte}")
    
    # Printing the final equation with each number as requested.
    print("\nFinal Equation for Maximal Size:")
    print(f"{total_bits_required} bits / {bits_in_a_byte} bits_per_byte = {size_in_bytes} bytes")
