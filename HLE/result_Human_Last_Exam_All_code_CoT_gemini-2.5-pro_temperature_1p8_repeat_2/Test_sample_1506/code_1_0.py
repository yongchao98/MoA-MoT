import sys

class PrimeDataStructure:
    """
    A data structure to efficiently check for primality and list primes
    up to a certain limit.
    """
    def __init__(self, limit=10000):
        """
        Initializes the data structure by pre-computing primes up to the limit
        using a Sieve of Eratosthenes.
        """
        if not isinstance(limit, int) or limit < 2:
            raise ValueError("Limit must be an integer greater than or equal to 2.")
        
        self.limit = limit
        # In a memory-optimized implementation, this would be a bit array.
        # For clarity and idiomatic Python, we use a list of booleans.
        self._is_prime = [True] * self.limit
        self._is_prime[0] = self._is_prime[1] = False

        # Sieve of Eratosthenes algorithm
        for i in range(2, int(self.limit**0.5) + 1):
            if self._is_prime[i]:
                # Mark all multiples of i as not prime
                for multiple in range(i * i, self.limit, i):
                    self._is_prime[multiple] = False

    def isprime(self, p):
        """
        Checks if a number p is prime in O(1) time.
        p must be less than the initialized limit.
        """
        if 0 <= p < self.limit:
            return self._is_prime[p]
        return False

    def primes(self, n):
        """
        Returns a list of all prime numbers less than or equal to n in O(n) time.
        n must be less than the initialized limit.
        """
        if n >= self.limit:
            n = self.limit - 1
            print(f"Warning: n is out of bounds, capping at {n}", file=sys.stderr)

        prime_list = []
        for i in range(n + 1):
            if self._is_prime[i]:
                prime_list.append(i)
        return prime_list

def main():
    """
    Main function to demonstrate the PrimeDataStructure and calculate its maximal size.
    """
    # The design covers numbers p < 10000 and n < 10000.
    # This requires storing primality information for numbers 0 through 9999.
    num_elements = 10000
    bits_per_byte = 8

    # The most memory-efficient data structure is a bit array, where each bit
    # represents whether a number is prime.
    max_size_in_bytes = num_elements / bits_per_byte
    
    print("Demonstrating the PrimeDataStructure:")
    checker = PrimeDataStructure(limit=num_elements)
    print(f"Is 9973 prime? {checker.isprime(9973)}")
    print(f"Is 100 prime? {checker.isprime(100)}")
    print(f"Primes up to 30: {checker.primes(30)}")
    print("-" * 20)

    print("Calculation for the maximal size of the data structure:")
    print(f"The most memory-efficient structure is a bit array to store primality for {num_elements} numbers (0-9999).")
    # This final print statement outputs each number in the equation as requested.
    print(f"Equation: {num_elements} bits / {bits_per_byte} bits per byte = {int(max_size_in_bytes)} bytes")

if __name__ == "__main__":
    main()
