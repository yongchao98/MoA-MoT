import math

class PrimeDataStructure:
    """
    A data structure to efficiently check for primality and list primes
    for numbers less than a given maximum value.
    """

    def __init__(self, max_val=10000):
        """
        Initializes the data structure by pre-computing primes using the
        Sieve of Eratosthenes.
        """
        self._max_val = max_val
        # Create a boolean list to store primality. Initialize all to True.
        self._sieve = [True] * self._max_val
        
        # 0 and 1 are not prime numbers
        if self._max_val > 0:
            self._sieve[0] = False
        if self._max_val > 1:
            self._sieve[1] = False
            
        # Mark multiples of primes as not prime
        for p in range(2, int(math.sqrt(self._max_val))):
            if self._sieve[p]:
                # Start marking from p*p
                for i in range(p * p, self._max_val, p):
                    self._sieve[i] = False

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be less than the initialized max_val.
        """
        if 0 <= p < self._max_val:
            return self._sieve[p]
        raise ValueError(f"Input p must be between 0 and {self._max_val - 1}")

    def primes(self, n):
        """
        Returns a list of all primes less than or equal to n in O(n) time.
        n must be less than the initialized max_val.
        """
        if 0 <= n < self._max_val:
            prime_list = []
            for p in range(2, n + 1):
                if self._sieve[p]:
                    prime_list.append(p)
            return prime_list
        raise ValueError(f"Input n must be between 0 and {self._max_val - 1}")

# --- Main Execution ---

# 1. Create an instance of the data structure
prime_checker = PrimeDataStructure(max_val=10000)

# 2. Demonstrate the operations
print("--- Demonstrating Operations ---")
p_to_check = 9973
is_p_prime = prime_checker.isprime(p_to_check)
print(f"Is {p_to_check} a prime? {is_p_prime}")

p_to_check = 9975
is_p_prime = prime_checker.isprime(p_to_check)
print(f"Is {p_to_check} a prime? {is_p_prime}")

n_for_list = 30
prime_list_up_to_n = prime_checker.primes(n_for_list)
print(f"Primes up to {n_for_list}: {prime_list_up_to_n}")
print("-" * 30)

# 3. Calculate and print the maximal size of the data structure
print("\n--- Maximal Size Calculation ---")
max_numbers = 10000
bits_per_number = 1
bits_in_a_byte = 8

total_bits = max_numbers * bits_per_number
total_bytes = total_bits // bits_in_a_byte

print("The most memory-efficient design uses a bit array.")
print(f"It needs to store primality for {max_numbers} numbers (0 to 9999).")
print(f"Each number's primality can be stored in {bits_per_number} bit.")
print(f"Total size in bytes = ({max_numbers} numbers * {bits_per_number} bit/number) / {bits_in_a_byte} bits/byte")
print(f"Final Equation: {total_bits} / {bits_in_a_byte} = {total_bytes} bytes")
