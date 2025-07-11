import random

def is_prime(n, k=7):
    """
    Miller-Rabin primality test.
    k is the number of rounds of testing to perform.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    
    # Write n as 2^r * d + 1
    d = n - 1
    r = 0
    while d % 2 == 0:
        r += 1
        d //= 2
    
    # Witness loop
    for _ in range(k):
        a = random.randrange(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def solve_puzzle():
    """
    Finds the 6th prime in the sequence from the digits of gamma (Euler-Mascheroni constant).
    """
    # The first 100 digits of the fractional part of the Euler-Mascheroni constant (gamma)
    gamma_digits = "577215664901532860606512090082402431042159045582327962148086411323381483135445210230324025534062"
    
    prime_count = 0
    found_primes = []
    
    for i in range(1, len(gamma_digits) + 1):
        num_str = gamma_digits[:i]
        num = int(num_str)
        
        if is_prime(num):
            prime_count += 1
            found_primes.append(num)
            if prime_count == 6:
                sixth_prime = num
                print(sixth_prime)
                return

solve_puzzle()