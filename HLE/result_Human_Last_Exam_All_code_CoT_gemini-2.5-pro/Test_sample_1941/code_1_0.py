import math

def is_prime(n):
    """
    Miller-Rabin primality test. This is a probabilistic test, but by using a
    specific set of bases, it can be made deterministic for numbers within a
    certain range, which is sufficient for this problem.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Using bases that are deterministic for n < 3,317,044,064,279,371
    # which is well beyond the expected size of numbers we'll encounter.
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    if n < 341550071728321:
        bases = [2, 3, 5, 7, 11, 13, 17]

    for a in bases:
        if a >= n:
            break
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def solve():
    """
    This function implements the plan to find the solution. It generates the
    required set of prime numbers and then counts how many of them will result
    in a Game of Life pattern that dies out.
    """
    
    # Step 1: Generate the first 1,000 prime numbers containing only '0' and '1' digits.
    # This set is referred to as C.
    primes_c = []
    limit = 1000
    i = 2  # Start checking from binary '10' (decimal 2), as '1' is not prime.
    
    while len(primes_c) < limit:
        # Generate a candidate number by interpreting the binary string of 'i' as a base-10 number.
        candidate_str = bin(i)[2:]
        candidate_num = int(candidate_str)

        if is_prime(candidate_num):
            primes_c.append(candidate_num)
        
        i += 1

    # Step 2: Count the members of C that will die out.
    # Based on our analysis, this corresponds to the numbers whose decimal representation
    # does not contain the substring "111".
    died_out_count = 0
    for p in primes_c:
        if '111' not in str(p):
            died_out_count += 1
    
    # Step 3: Print the final count as the result.
    # The instruction "output each number in the final equation" is interpreted as
    # printing the single, final numerical answer of this calculation.
    print(died_out_count)

solve()