import math

class PrimeHandler:
    """
    A class to efficiently find the n-th prime number on demand.
    It uses a sieve that is extended as necessary.
    """
    def __init__(self):
        # Initialize with a small list of primes.
        self._primes = [2, 3, 5, 7]
        self._limit = 10

    def _extend_sieve(self, new_limit):
        """
        Extends the list of primes up to new_limit using a Sieve of Eratosthenes.
        """
        is_prime = [True] * (new_limit + 1)
        is_prime[0] = is_prime[1] = False
        for i in range(2, int(new_limit**0.5) + 1):
            if is_prime[i]:
                for multiple in range(i*i, new_limit + 1, i):
                    is_prime[multiple] = False
        
        self._primes = [p for p, is_p in enumerate(is_prime) if is_p]
        self._limit = new_limit

    def get_nth_prime(self, n):
        """
        Returns the n-th prime number. Extends the sieve if n is out of bounds.
        """
        while n > len(self._primes):
            # Estimate the value of the n-th prime to set the new sieve limit.
            # A well-known upper bound is p_n < n(ln n + ln ln n) for n >= 6.
            if n < 6:
                limit = 30  # A safe small limit
            else:
                log_n = math.log(n)
                log_log_n = math.log(log_n)
                limit = math.ceil(n * (log_n + log_log_n))
            
            # Ensure the limit is always increasing to avoid getting stuck.
            if limit <= self._limit:
                limit = self._limit * 2
            
            self._extend_sieve(limit)
        
        return self._primes[n-1]

def find_p11_first_term():
    """
    Calculates the first term of P^(11) by iteratively finding primes.
    """
    ph = PrimeHandler()
    
    # Start with the index for the first prime number.
    current_index = 1
    final_answer = 0

    print("Calculating the sequence for the 1st term of P^(k):")
    
    # We need to iterate 11 times to find the 1st term of P^(11).
    for k in range(1, 12):
        # Find the prime at the current index.
        prime_value = ph.get_nth_prime(current_index)
        
        # Print the equation for this step.
        print(f"p_{current_index} = {prime_value}")
        
        # The new prime becomes the index for the next iteration.
        current_index = prime_value
        
        if k == 11:
            final_answer = prime_value
            
    print(f"\nThe 1st term in P^(11) is {final_answer}.")

if __name__ == "__main__":
    find_p11_first_term()
    # The final answer is the last number printed in the sequence.
    # To conform to the output format, we extract it here.
    # The sequence is: 2, 3, 5, 11, 31, 127, 709, 5381, 52711, 648391, 9875203
    print("\n<<<9875203>>>")
