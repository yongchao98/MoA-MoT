def solve():
    """
    Calculates the global labeling number of the graph K_1,100.
    The number is determined to be 2^99.
    """
    base = 2
    exponent = 100 - 1
    
    # Calculate the result
    result = base ** exponent
    
    # Print the equation and the result
    print(f"The global labeling number of K_1,100 is given by the equation:")
    print(f"{base} ** {exponent} = {result}")

solve()