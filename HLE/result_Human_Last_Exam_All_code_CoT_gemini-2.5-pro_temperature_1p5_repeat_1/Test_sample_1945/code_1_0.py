import math

def calculate_limit():
    """
    This function calculates the exact value of the limit lim_{n->inf} Pr(n).
    The formula for the limit is 3 / (2 * pi^2).
    """
    numerator = 3
    denominator_part1 = 2
    denominator_part2 = math.pi ** 2
    
    limit_value = numerator / (denominator_part1 * denominator_part2)
    
    # The instruction requires printing the numbers in the final equation.
    print(f"The equation for the limit is: {numerator} / ({denominator_part1} * {denominator_part2})")
    print("The exact value of the limit is:")
    print(limit_value)

calculate_limit()