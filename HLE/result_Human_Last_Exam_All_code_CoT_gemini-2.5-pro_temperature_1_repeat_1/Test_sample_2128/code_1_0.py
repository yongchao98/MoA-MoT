import math

def solve_sequence_problem():
    """
    Calculates the value of 1/p_1000 based on the derived formula.
    The formula is 1/p_n = 4 * cos^2(pi / (n + 2)).
    For n = 1000, this becomes 1/p_1000 = 4 * cos^2(pi / 1002).
    """
    n = 1000
    
    # The numbers in the final equation: 1/p_1000 = 4 * cos^2(pi / 1002)
    constant_factor = 4
    denominator_in_arg = n + 2
    
    # Calculate the value
    angle = math.pi / denominator_in_arg
    cos_val = math.cos(angle)
    result = constant_factor * (cos_val ** 2)
    
    # Print the equation and the numbers involved
    print(f"The final equation is: 1/p_{n} = {constant_factor} * cos^2(pi / {denominator_in_arg})")
    
    # Print the final calculated value
    print(f"The value of 1/p_{n} for n = {n} is: {result}")

solve_sequence_problem()
