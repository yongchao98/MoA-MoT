import math

def solve_pm():
    """
    This function calculates the probability P_m based on the logic derived.
    You can change the value of m to explore different cases.
    """
    # You can change the value of m here.
    m = 4 

    print(f"Analyzing the problem for m = {m}")

    # The denominator of the probability is the total number of ways to choose
    # 2 items from a set of 4m+2, which is C(4*m + 2, 2).
    # C(n, 2) = n * (n - 1) / 2
    # C(4*m + 2, 2) = ((4*m + 2) * (4*m + 1)) / 2 = (2*m + 1) * (4*m + 1)
    
    total_pairs_denominator_expression = "(2*m + 1)*(4*m + 1)"
    
    if m % 2 == 1:
        # Case when m is odd
        numerator = 3
        print(f"Since m is odd, the number of valid (i, j) pairs that allow the remaining items to be partitioned is {numerator}.")
        final_equation = f"P_m = {numerator} / {total_pairs_denominator_expression}"
        print(f"The probability formula is: {final_equation}")
        
    else:
        # Case when m is even
        numerator = 4
        print(f"Since m is even, the number of valid (i, j) pairs is {numerator}.")
        final_equation = f"P_m = {numerator} / {total_pairs_denominator_expression}"
        print(f"The probability formula is: {final_equation}")

    # Now, calculate the numerical value for the given m
    val_2 = 2
    val_1_1 = 1
    val_4 = 4
    val_1_2 = 1
    
    denominator = (val_2 * m + val_1_1) * (val_4 * m + val_1_2)
    
    # Simplify the fraction
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor

    print(f"\nFor the specific value m = {m}:")
    print(f"Total pairs = {denominator}")
    print(f"The probability P_{m} = {numerator}/{denominator} = {simplified_num}/{simplified_den}")
    
    print("\nThe numbers that form the general final equation are:")
    print(f"Numerator: {numerator}")
    print(f"Numbers in the denominator expression (2*m+1)*(4*m+1): {val_2}, {val_1_1}, {val_4}, {val_1_2}")

solve_pm()