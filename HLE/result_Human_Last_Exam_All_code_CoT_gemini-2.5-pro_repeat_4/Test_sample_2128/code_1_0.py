import math

def solve():
    """
    This function calculates the value of 1/p_1000 based on the derived formula.
    """
    n = 1000
    
    # The formula for 1/p_n is 4 * cos^2(pi / (n + 2)).
    # We substitute n = 1000 into the formula.
    
    numerator = 4
    n_plus_2 = n + 2
    
    # Calculate the value
    result = numerator * (math.cos(math.pi / n_plus_2))**2
    
    # Print the final equation with all the numbers and the result
    print(f"The value of 1/p_{n} is given by the equation:")
    print(f"1/p_{n} = 4 * cos^2(pi / (n + 2))")
    print(f"For n = {n}:")
    print(f"1/p_{n} = 4 * cos^2(pi / {n_plus_2})")
    print(f"The calculated value is: {result}")

solve()