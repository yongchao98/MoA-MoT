import random

def is_prime(n, k=5):
    """
    Tests if a number is prime using the Miller-Rabin primality test.
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
        d //= 2
        r += 1

    # Witness loop for k rounds
    for _ in range(k):
        a = random.randrange(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        
        is_composite = True
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                is_composite = False
                break
        if is_composite:
            return False
            
    return True

def solve_conway_primes():
    """
    Finds the first 1,000 primes with only digits 0 and 1,
    and counts how many would die out in Conway's Game of Life.
    """
    limit = 1000
    primes_found_count = 0
    dying_patterns_count = 0
    
    # Start with integer 2, whose binary representation is '10'
    i = 2
    
    print("Finding the first 1,000 primes containing only digits 0 and 1...")
    print("Then, counting how many of them die out in Conway's Game of Life.")
    print("(A pattern dies out if its number does not contain the substring '111')")
    
    while primes_found_count < limit:
        # Generate the number string (e.g., for i=5, binary_str is "101")
        binary_str = bin(i)[2:]
        
        # Convert the string to an integer (e.g., "101" -> 101)
        num_to_test = int(binary_str)
        
        if is_prime(num_to_test):
            primes_found_count += 1
            
            # The pattern dies out if and only if it does NOT contain '111'
            if "111" not in binary_str:
                dying_patterns_count += 1
        
        i += 1
        
    print("\n--- Calculation Complete ---")
    print("Out of the first 1,000 prime numbers containing only digits 0 and 1:")
    # Per instructions, print the final number as the "result of the equation"
    print(f"The number of members that will ultimately die out is: {dying_patterns_count}")

solve_conway_primes()