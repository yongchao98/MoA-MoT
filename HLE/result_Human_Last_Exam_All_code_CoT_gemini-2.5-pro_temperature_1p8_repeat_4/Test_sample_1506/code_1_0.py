import math

class PrimeDataStructure:
    """
    A memory-efficient data structure for prime number operations up to a certain limit.
    It uses a bit array storing primality only for odd numbers to save space.
    """
    def __init__(self, limit=9999):
        """
        Initializes the data structure by sieving for primes up to the limit.
        """
        self.limit = limit
        # We only store primality for odd numbers. An odd number p is mapped to index (p-1)//2.
        # e.g., 1 -> 0, 3 -> 1, 5 -> 2, ..., 9999 -> 4999.
        # The number of bits needed is (limit + 1) / 2.
        num_bits = (self.limit + 1) // 2
        self.num_bytes = math.ceil(num_bits / 8)
        
        # Initialize a bytearray with all bits set to 1 (True).
        self._sieve = bytearray([0xFF] * self.num_bytes)

        # 1 is not prime. Its index is (1-1)//2 = 0.
        self._set_bit(0, False)

        # Perform the Sieve of Eratosthenes for odd numbers.
        # Start with the first odd prime, 3.
        # Iterate up to the square root of the limit.
        for p in range(3, int(self.limit**0.5) + 1, 2):
            # If p is prime (its bit is still 1)...
            if self._get_bit((p - 1) // 2):
                # ...then sieve out its multiples.
                # Start sieving from p*p.
                # Multiples are p*p, p*p + 2p, p*p + 4p, etc. (all are odd).
                for multiple in range(p * p, self.limit + 1, 2 * p):
                    self._set_bit((multiple - 1) // 2, False)

    def _get_bit(self, index):
        """Gets the bit value at a given index in the sieve."""
        byte_index = index // 8
        bit_index = index % 8
        return (self._sieve[byte_index] >> bit_index) & 1

    def _set_bit(self, index, value):
        """Sets the bit value at a given index in the sieve."""
        byte_index = index // 8
        bit_index = index % 8
        if value:
            # Set bit to 1
            self._sieve[byte_index] |= (1 << bit_index)
        else:
            # Set bit to 0
            self._sieve[byte_index] &= ~(1 << bit_index)

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be less than 10000.
        """
        if not 0 <= p <= self.limit:
            raise ValueError(f"Input {p} is out of the supported range [0, {self.limit}]")
        
        if p == 2:
            return True
        if p < 2 or p % 2 == 0:
            return False
            
        # For odd p, look up its primality in the sieve.
        return bool(self._get_bit((p - 1) // 2))

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n in O(n) time.
        n must be less than 10000.
        """
        if not 0 <= n <= self.limit:
            raise ValueError(f"Input {n} is out of the supported range [0, {self.limit}]")
        
        prime_list = []
        if n >= 2:
            prime_list.append(2)
        
        # Iterate through odd numbers up to n and check for primality.
        for p in range(3, n + 1, 2):
            if self.isprime(p):
                prime_list.append(p)
        return prime_list

def main():
    """Main function to demonstrate the data structure and print its size."""
    # The problem specifies checks for p < 10000 and n < 10000,
    # so the highest number we need to consider is 9999.
    max_number = 9999
    
    # Initialize the prime data structure.
    prime_db = PrimeDataStructure(limit=max_number)

    # --- Demonstrations ---
    print(f"Is 9973 prime? {prime_db.isprime(9973)}")
    print(f"Is 9999 prime? {prime_db.isprime(9999)}")
    print(f"Primes up to 30: {prime_db.primes(30)}")
    
    # --- Maximal Size Calculation ---
    # The maximal size is determined by the bit array holding the primality data.
    # The size is fixed on initialization.
    num_supported = max_number + 1
    num_odd_numbers = num_supported // 2
    bits_per_byte = 8
    
    # The bytearray stores the data. We print the equation for its size calculation.
    print("\n--- Maximal Size Calculation ---")
    print(f"The size of the core data structure is determined by the number of bits needed.")
    print(f"Numbers to store primality for: Odd numbers from 1 to {max_number}")
    print(f"Total odd numbers = {num_odd_numbers}")
    print(f"Total bits needed = {num_odd_numbers}")
    print(f"Equation for size in bytes: ceil( Total bits / Bits per byte )")
    print(f"Size = ceil( {num_odd_numbers} / {bits_per_byte} ) = {prime_db.num_bytes} bytes.")
    print("(Note: This is the size of the raw data. The Python object itself has additional overhead.)")

if __name__ == "__main__":
    main()
