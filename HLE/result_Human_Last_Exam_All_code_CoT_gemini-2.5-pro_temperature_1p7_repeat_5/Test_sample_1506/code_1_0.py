import math

class PrimeDataStructure:
    """
    A data structure to efficiently check for primality and list primes
    up to a certain limit.
    """
    def __init__(self, max_num=10000):
        """
        Initializes the data structure by building a prime sieve.
        """
        self._max_num = max_num
        # Create a boolean list to store primality.
        # Initially, assume all numbers are prime.
        self._sieve = [True] * self._max_num
        
        # 0 and 1 are not prime numbers.
        if self._max_num > 0:
            self._sieve[0] = False
        if self._max_num > 1:
            self._sieve[1] = False
        
        # Use Sieve of Eratosthenes to mark non-prime numbers.
        # We only need to check up to the square root of max_num.
        for number in range(2, int(math.sqrt(self._max_num))):
            if self._sieve[number]:
                # Mark all multiples of the current prime number as not prime.
                # Start marking from number*number.
                for multiple in range(number * number, self._max_num, number):
                    self._sieve[multiple] = False

    def isprime(self, p):
        """
        Checks if a number p is prime.
        Time Complexity: O(1)
        """
        if not 0 <= p < self._max_num:
            raise ValueError(f"Input {p} is out of the supported range [0, {self._max_num-1}]")
        return self._sieve[p]

    def primes(self, n):
        """
        Returns a list of all primes p <= n.
        Time Complexity: O(n)
        """
        if not 0 <= n < self._max_num:
            raise ValueError(f"Input {n} is out of the supported range [0, {self._max_num-1}]")
        
        prime_list = []
        # Iterate from 2 to n and collect all numbers marked as prime.
        for number in range(2, n + 1):
            if self._sieve[number]:
                prime_list.append(number)
        return prime_list

# --- Demonstration ---

# 1. Create an instance of the data structure.
#    This pre-computes all primes up to 9999.
prime_checker = PrimeDataStructure(10000)

# 2. Use the `isprime` method (O(1) time).
p1 = 29
p2 = 30
print(f"Is {p1} prime? {prime_checker.isprime(p1)}")
print(f"Is {p2} prime? {prime_checker.isprime(p2)}")

# 3. Use the `primes` method (O(n) time).
n = 30
prime_list_up_to_30 = prime_checker.primes(n)
print(f"Primes up to {n}: {prime_list_up_to_30}")

# 4. Calculate and print the maximal size of the data structure.
# The structure stores a boolean for 10,000 numbers (0 to 9999).
# Most efficient storage uses 1 bit per number.
total_numbers = 10000
bits_per_byte = 8
max_size_in_bytes = total_numbers / bits_per_byte

# The problem asks for the final answer in the format <<<answer>>>
# Here we print the equation for clarity before printing the final result.
print("\nCalculating the maximal size of the data structure in bytes:")
print(f"The data structure needs to store a boolean for each of the {total_numbers} numbers.")
print(f"The most memory-efficient way is to use 1 bit per number.")
print(f"Total bits = {total_numbers}")
print(f"Size in bytes = Total bits / bits per byte")
print(f"Size in bytes = {total_numbers} / {bits_per_byte} = {int(max_size_in_bytes)}")
<<<1250>>>