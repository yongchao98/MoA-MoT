import math

def solve():
    """
    Calculates the value of 1/p_1000 based on the derived formula.
    """
    n = 1000
    
    # The formula for 1/p_n is 4 * cos^2(pi / (n + 2))
    four = 4
    n_plus_2 = n + 2
    
    # Calculate the value
    result = four * (math.cos(math.pi / n_plus_2))**2
    
    # Print the equation with the numbers plugged in, and the final result.
    print(f"The formula for 1/p_n is: 4 * cos^2(pi / (n + 2))")
    print(f"For n = {n}, the equation is:")
    print(f"1/p_{n} = {four} * cos^2(pi / {n_plus_2})")
    print(f"The calculated value is: {result}")

solve()