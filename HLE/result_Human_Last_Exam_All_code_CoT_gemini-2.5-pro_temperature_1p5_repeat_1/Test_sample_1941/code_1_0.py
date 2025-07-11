import math

def is_prime(n):
    """
    Efficiently checks if a number is prime.
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

def solve_game_of_life_primes():
    """
    Finds the first 1,000 prime numbers composed only of digits 0 and 1,
    and counts how many of them would "die out" in Conway's Game of Life.
    
    A one-cell-thick pattern dies out if and only if its decimal representation
    does not contain the substring '111'.
    """
    target_count = 1000
    found_primes = 0
    dies_out_count = 0
    
    # We generate numbers with only digits 0 and 1 by using the binary 
    # representation of an incrementing integer 'i'.
    i = 2  # Start at 2 because bin(1)='1' -> n=1 is not prime.
    while found_primes < target_count:
        # Step 1: Generate a candidate number string (e.g., for i=5, s='101')
        s = bin(i)[2:]
        
        # Convert the string representation to an integer
        n = int(s)
        
        # Step 2: Check if the number is prime
        if is_prime(n):
            # We found a member of Set C
            found_primes += 1
            
            # Step 3: Check the Game of Life condition for dying out
            # The pattern dies if it does NOT contain a '111' block.
            if '111' not in s:
                dies_out_count += 1
        
        i += 1
        
    print(f"Out of the first {target_count} primes in set C, the number that will ultimately die out is:")
    print(dies_out_count)

solve_game_of_life_primes()
<<<819>>>