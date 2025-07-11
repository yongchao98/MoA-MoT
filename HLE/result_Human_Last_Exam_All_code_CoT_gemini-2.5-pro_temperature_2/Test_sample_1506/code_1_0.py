import math

class PrimeSieve:
    """
    A memory-efficient data structure for primality testing and generation
    using a bit-packed sieve of odds. It provides O(1) for isprime(p)
    and O(n) for primes(n).
    """

    def __init__(self, limit):
        """
        Initializes the sieve for all numbers less than limit.
        """
        if limit < 3:
            raise ValueError("Limit must be at least 3.")
        self.limit = limit
        
        # We only need to sieve odd numbers >= 3. The index `i` in the sieve
        # corresponds to the odd number `2*i + 3`.
        # The number of bits needed is for odds from 3 up to limit-1.
        sieve_size_in_bits = (limit - 1 - 3) // 2 + 1 if limit > 3 else 0
        
        num_bytes = (sieve_size_in_bits + 7) // 8
        self._sieve_data = bytearray([0xff] * num_bytes)  # 1 represents prime

        # --- Helper functions to manipulate bits based on the number ---
        def _clear_bit(p):
            # p must be an odd number >= 3
            idx = (p - 3) // 2
            self._sieve_data[idx >> 3] &= ~(1 << (idx & 7))

        def _test_bit(p):
            # p must be an odd number >= 3
            idx = (p - 3) // 2
            return (self._sieve_data[idx >> 3] >> (idx & 7)) & 1

        # --- Perform Sieve of Eratosthenes (for odds) ---
        # We only need to iterate up to sqrt(limit)
        for p in range(3, int(limit**0.5) + 1, 2):
            if _test_bit(p):
                # p is a prime. Mark all its odd multiples as not prime.
                # Start from p*p, as smaller multiples are already handled.
                # The step is 2*p to only hit odd multiples (p*p, p*(p+2), etc.).
                for multiple in range(p * p, limit, 2 * p):
                    _clear_bit(multiple)

    def isprime(self, p):
        """
        Checks if a number p is prime. Time complexity: O(1).
        Assumes p is within the initialized limit.
        """
        if p < 2 or p >= self.limit:
            return False
        if p == 2:
            return True
        if p % 2 == 0:
            return False
        
        # Now p is an odd number >= 3. Look it up in our sieve.
        idx = (p - 3) // 2
        return (self._sieve_data[idx >> 3] >> (idx & 7)) & 1 == 1

    def primes(self, n):
        """
        Returns a list of all prime numbers p <= n. Time complexity: O(n).
        Assumes n is within the initialized limit.
        """
        if n >= self.limit:
            n = self.limit - 1
            
        if n < 2:
            return []
        
        prime_list = [2]  # Start with the only even prime
        
        # Iterate through odd numbers and check for primality
        for p in range(3, n + 1, 2):
            if self.isprime(p):
                prime_list.append(p)
                
        return prime_list


if __name__ == '__main__':
    LIMIT = 10000

    # You can instantiate and use the class like this:
    # sieve = PrimeSieve(LIMIT)
    # print(f"Is 9999 prime? {sieve.isprime(9999)}")
    # print(f"Is 9991 prime? {sieve.isprime(9991)}")
    # print(f"Primes up to 30: {sieve.primes(30)}")

    print("This script designs a data structure and calculates its maximal size.")
    print("The design is a memory-optimized 'sieve of odds' stored in a bit array.\n")

    print("--- Calculation of the Maximal Data Structure Size ---")

    max_num_to_check = LIMIT - 1
    
    print("\n1. Determine the number of bits required:")
    print(f"The structure handles primes up to p < {LIMIT}, i.e., up to {max_num_to_check}.")
    print("We only need to store primality for odd numbers from 3 upwards.")
    print("Number of bits = (max_num - 3) / 2 + 1")
    bits_eq_step1 = f"Number of bits = ({max_num_to_check} - 3) / 2 + 1"
    print(bits_eq_step1)
    num_of_odds = (max_num_to_check - 3) // 2 + 1
    bits_eq_step2 = f"Number of bits = {max_num_to_check - 3} / 2 + 1 = {num_of_odds}"
    print(bits_eq_step2)
    
    print("\n2. Convert bits to bytes:")
    print("Size in bytes = ceil(Number of bits / 8)")
    bytes_eq_step1 = f"Size in bytes = ceil({num_of_odds} / 8)"
    print(bytes_eq_step1)
    bytes_float = num_of_odds / 8
    bytes_eq_step2 = f"Size in bytes = ceil({bytes_float})"
    print(bytes_eq_step2)
    final_size_bytes = math.ceil(bytes_float)
    bytes_eq_step3 = f"Size in bytes = {final_size_bytes}"
    print(bytes_eq_step3)

    print(f"\nThe maximal size of this data structure is {final_size_bytes} bytes.")
    print("<<<625>>>")