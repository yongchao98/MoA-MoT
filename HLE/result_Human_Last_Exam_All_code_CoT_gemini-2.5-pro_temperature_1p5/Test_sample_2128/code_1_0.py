import math

def solve():
    """
    Calculates the value of 1/p_1000 based on the derived formula.
    """
    n = 1000
    
    # The final expression for 1/p_n is 4 * cos^2(pi / (n + 2)).
    # We substitute n = 1000 into this expression.
    
    # The number n in the problem statement.
    n_val = 1000
    
    # The number 2 from the formula n + 2.
    offset = 2
    
    # The constant factor in the expression.
    factor = 4
    
    # Calculate the denominator inside the cosine function.
    denominator = n_val + offset
    
    # The argument of the cosine function.
    angle = math.pi / denominator
    
    # Calculate the final value.
    result = factor * (math.cos(angle) ** 2)

    # Output the steps as requested, showing each number in the final equation.
    print(f"The expression for 1/p_n is: 1/p_n = {factor} * cos^2(pi / (n + {offset}))")
    print(f"For n = {n_val}, the expression is: 1/p_{n_val} = {factor} * cos^2(pi / {denominator})")
    print(f"The calculated value is: {result}")

solve()