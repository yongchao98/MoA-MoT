import random

def power(a, b, m):
    """
    Calculates (a^b) % m using modular exponentiation.
    """
    res = 1
    a %= m
    while b > 0:
        if b % 2 == 1:
            res = (res * a) % m
        a = (a * a) % m
        b //= 2
    return res

def is_prime(n, k=10):
    """
    Miller-Rabin primality test. Returns True if n is probably prime.
    k is the number of rounds of testing.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False

    # Write n-1 as 2^r * d
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    
    # Test for k random bases
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = power(a, d, n)
        
        if x == 1 or x == n - 1:
            continue
        
        is_composite = True
        for _ in range(r - 1):
            x = power(x, 2, n)
            if x == n - 1:
                is_composite = False
                break
        
        if is_composite:
            return False
            
    return True # Probably prime

def solve_game_of_life_primes():
    """
    Finds the first 1,000 primes with only 0/1 digits and counts
    how many would die out in Conway's Game of Life.
    """
    primes_c = []
    dying_count = 0
    num_generator = 2 # Start with i=2 to generate candidate "10"

    while len(primes_c) < 1000:
        # Generate a candidate number n whose decimal digits are only 0s and 1s.
        # This is done by interpreting the binary representation of an integer as a decimal number.
        # e.g., i=3 -> bin(3)="11" -> n=11. i=5 -> bin(5)="101" -> n=101.
        s = bin(num_generator)[2:]
        
        # Optimization 1: A number with 0/1 digits ending in 0 is divisible by 10, so not prime.
        if s.endswith('0'):
            num_generator += 1
            continue
            
        # Optimization 2: A number is divisible by 3 if its digit sum is divisible by 3.
        # Here, the digit sum is just the count of '1's.
        # Since all our numbers are > 3, this check is safe.
        if s.count('1') % 3 == 0:
            num_generator += 1
            continue

        n = int(s)
        
        # Test if the generated number is prime.
        if is_prime(n):
            primes_c.append(n)
            
            # A 1D Game of Life pattern dies out if its string representation does not contain "111".
            if "111" not in s:
                dying_count += 1

        num_generator += 1

    surviving_count = 1000 - dying_count
    
    # The final equation shows the total set size, the non-dying count, and the dying count.
    print(f"{len(primes_c)} - {surviving_count} = {dying_count}")

solve_game_of_life_primes()