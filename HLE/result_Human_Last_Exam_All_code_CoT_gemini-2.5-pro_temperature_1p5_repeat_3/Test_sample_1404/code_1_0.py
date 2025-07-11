import math

def calculate_arrangements():
    """
    Calculates the number of arrangements for both large and small packages
    and determines if the small package has 1260 times more arrangements.
    """
    # Step 1: Calculate arrangements for the large package
    # 9 distinct red candles for 9 horizontal positions -> 9!
    # 9 distinct green candles for 9 vertical positions -> 9!
    # Total arrangements for large package = 9! * 9!
    factorial_9 = math.factorial(9)
    arrangements_large = factorial_9 * factorial_9

    print(f"Number of arrangements for the large package:")
    print(f"This is the permutation of 9 red candles multiplied by the permutation of 9 green candles.")
    print(f"Equation: 9! * 9! = {factorial_9} * {factorial_9} = {arrangements_large}")
    print("-" * 30)

    # Step 2: Calculate arrangements for the small package
    # 16 distinct candles for 16 horizontal positions -> 16!
    arrangements_small = math.factorial(16)
    
    print(f"Number of arrangements for the small package:")
    print(f"This is the permutation of 16 distinct candles.")
    print(f"Equation: 16! = {arrangements_small}")
    print("-" * 30)

    # Step 3: Calculate the ratio
    if arrangements_large > 0:
        ratio = arrangements_small / arrangements_large
    else:
        ratio = 0
        
    print(f"Ratio of small package arrangements to large package arrangements:")
    print(f"Equation: A_small / A_large = {arrangements_small} / {arrangements_large}")
    print(f"Calculated Ratio = {ratio:.4f}")
    print("-" * 30)
    
    # Step 4: Compare and conclude
    target_ratio = 1260
    is_true = math.isclose(ratio, target_ratio)
    
    print(f"Is it true that the number of arrangements for the small packages is {target_ratio} times greater than for the large packages?")
    print(f"The calculated ratio is approximately {ratio:.4f}, which is not equal to {target_ratio}.")
    print(f"Therefore, the statement is {is_true}.")

calculate_arrangements()