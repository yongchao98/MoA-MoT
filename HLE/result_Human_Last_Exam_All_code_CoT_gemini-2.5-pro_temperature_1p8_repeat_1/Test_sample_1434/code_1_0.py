import math

def calculate_sum():
    """
    Calculates the value of the sum described in the problem.
    The sum is broken down into two parts.
    """
    
    # Part 1: 2 * sum_{n=1 to inf} (1/4)^n
    # This is a geometric series: 2 * ( (1/4) / (1 - 1/4) )
    part1_val = 2 * ( (1/4) / (1 - 1/4) )

    # Part 2: 2 * sum_{n=1 to inf} 2^{v_2(n)} / 4^n
    # This sum, using a known identity, evaluates to 2 * ( (1/4) / (1-1/4) )
    # Let T(x) = sum_{n=1 to inf} 2^{v_2(n)} * x^n
    # We need 2 * T(1/4).
    # The identity is 2*T(x) = 2*x / (1-x)
    x = 1/4
    part2_val = 2 * x / (1 - x)

    # The total sum is Part 1 + Part 2
    total_sum = part1_val + part2_val
    
    print(f"The first part of the sum is 2 * Sum(1/4^n) = {part1_val}")
    print(f"The second part of the sum is 2 * Sum(2^v_2(n)/4^n) = {part2_val}")
    print(f"The total sum is {part1_val} + {part2_val} = {total_sum}")

calculate_sum()