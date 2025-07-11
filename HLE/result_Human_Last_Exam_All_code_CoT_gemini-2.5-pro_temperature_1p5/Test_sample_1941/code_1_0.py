import collections
import math

def is_prime(n):
    """
    A primality test function optimized for large numbers.
    It uses the 6k +/- 1 optimization.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    d = 5
    while d * d <= n:
        if n % d == 0 or n % (d + 2) == 0:
            return False
        d += 6
    return True

def solve():
    """
    This function finds the first 1,000 prime numbers containing only digits 0 and 1,
    and then determines how many of them would die out in Conway's Game of Life.
    """
    # C will store the first 1,000 binary-digit prime strings.
    C = []
    # A queue to generate binary-digit numbers in a breadth-first search manner.
    q = collections.deque(['1'])

    # Step 1: Find the first 1,000 binary-digit primes.
    while len(C) < 1000:
        s = q.popleft()

        # Enqueue the next potential candidates.
        q.append(s + '0')
        q.append(s + '1')
        
        # Optimization: A prime number (except 2) cannot be even.
        # So we only need to test numbers ending in '1'. Also, '1' is not prime.
        if s[-1] == '0' or s == '1':
            continue

        # Optimization: A number is divisible by 3 if the sum of its digits is.
        # For binary-digit numbers, this sum is just the count of '1's.
        if s.count('1') % 3 == 0:
            continue
            
        n = int(s)
        if is_prime(n):
            C.append(s)

    # Step 2: Count members of C that will die out.
    # This happens if their string representation does not contain "111".
    die_out_count = 0
    for prime_str in C:
        if "111" not in prime_str:
            die_out_count += 1
            
    survive_count = 1000 - die_out_count

    # The problem asks for the final answer in an equation format,
    # outputting each number involved. We present the total number of primes,
    # the number that survive, and the resulting number that die out.
    print(f"1000 - {survive_count} = {die_out_count}")

solve()
<<<838>>>