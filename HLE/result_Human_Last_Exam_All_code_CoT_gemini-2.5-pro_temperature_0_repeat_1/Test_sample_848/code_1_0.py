import math

def solve():
    """
    This function calculates the integer part of the limit expression derived from the problem.
    The final expression for the limit is L = 2 / ln(phi), where phi is the golden ratio-like constant
    from the recurrence relation of the solutions.
    The value to be computed is floor(10^4 * L).
    """

    # The final equation is 10000 * (2 / ln((13 + sqrt(165)) / 2))
    # The numbers in the equation are:
    val_10000 = 10000
    val_2_numerator = 2
    val_13 = 13
    val_165 = 165
    val_2_denominator = 2
    
    # We print the numbers as requested by the prompt.
    print(f"The numbers in the final equation are:")
    print(f"Multiplier: {val_10000}")
    print(f"Numerator: {val_2_numerator}")
    print(f"Term in ln: {val_13}")
    print(f"Term under square root: {val_165}")
    print(f"Denominator in ln: {val_2_denominator}")
    
    # Calculate phi
    phi = (val_13 + math.sqrt(val_165)) / val_2_denominator
    
    # Calculate the limit value L
    limit_L = val_2_numerator / math.log(phi)
    
    # Calculate the final result
    result = val_10000 * limit_L
    
    # Print the integer part of the result
    print("\nThe integer part of the final result is:")
    print(int(result))

solve()