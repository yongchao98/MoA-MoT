import math
from collections import deque

def is_prime(n: int) -> bool:
    """Checks if a number is prime."""
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    # Check for odd divisors from 3 up to sqrt(n).
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True

def solve():
    """
    Finds the first 1,000 primes with only 0/1 digits,
    counts how many will die out in Conway's Game of Life,
    and prints the result.
    """
    # Find the first 1,000 primes with only 0/1 digits.
    set_C = []
    # Use a queue to generate numbers with only 0/1 digits in ascending order.
    q = deque(['1']) 

    while len(set_C) < 1000:
        current_str = q.popleft()
        
        # The number 1 is not prime. Primes with only 0/1 digits start with 11.
        if len(current_str) > 1:
            current_num = int(current_str)
            if is_prime(current_num):
                set_C.append(current_num)
        
        # Add next candidates to the queue for breadth-first generation.
        q.append(current_str + '0')
        q.append(current_str + '1')

    # Count how many members of C will die out.
    # A pattern dies out if its string representation does not contain "111".
    die_out_count = 0
    for p in set_C:
        if "111" not in str(p):
            die_out_count += 1

    # As requested, we "output each number in the final equation" by representing
    # each prime that dies out as a '1' in a summation.
    if die_out_count > 0:
        equation_parts = ["1"] * die_out_count
        equation_str = " + ".join(equation_parts)
        print(f"{equation_str} = {die_out_count}")
    else:
        # This case is unlikely but included for completeness.
        print("0 = 0")

solve()