import math

def solve_candle_arrangement():
    """
    Calculates and compares the number of arrangements for two types of candle packages.
    """
    
    # Numbers for the large package
    num_red_large = 9
    num_green_large = 9
    
    # Numbers for the small package
    num_total_small = 16
    
    # Factor for comparison
    comparison_factor = 1260
    
    # --- Step 1: Calculate arrangements for the large package ---
    print("Step 1: Calculating arrangements for the large package (N_large)")
    print(f"There are {num_red_large} distinct red candles for {num_red_large} horizontal positions.")
    print(f"There are {num_green_large} distinct green candles for {num_green_large} vertical positions.")
    
    # Calculate 9!
    n9_factorial = math.factorial(num_red_large)
    print(f"The number of permutations for each color is {num_red_large}! = {n9_factorial}")
    
    # Calculate N_large = 9! * 9!
    n_large = n9_factorial * n9_factorial
    print(f"Total arrangements for the large package = {num_red_large}! * {num_green_large}! = {n9_factorial} * {n9_factorial} = {n_large}\n")

    # --- Step 2: Calculate arrangements for the small package ---
    print("Step 2: Calculating arrangements for the small package (N_small)")
    print(f"There are {num_total_small} distinct candles to be placed in {num_total_small} horizontal positions.")
    
    # Calculate 16!
    n16_factorial = math.factorial(num_total_small)
    print(f"The number of permutations is {num_total_small}! = {n16_factorial}")
    
    n_small = n16_factorial
    print(f"Total arrangements for the small package = {n16_factorial}\n")

    # --- Step 3: Compare the arrangements ---
    print("Step 3: Calculating the ratio (N_small / N_large)")
    print(f"The ratio is {num_total_small}! / ({num_red_large}! * {num_green_large}!)")
    
    # Calculate the ratio
    ratio = n_small / n_large
    print(f"Ratio = {n_small} / {n_large} = {ratio}\n")

    # --- Step 4: Final Conclusion ---
    print("Step 4: Checking if the ratio is equal to 1260")
    is_true = (ratio == comparison_factor)
    print(f"The calculated ratio is {ratio}.")
    print(f"The proposed ratio is {comparison_factor}.")
    print(f"Is it true that the number of arrangements for the small packages is {comparison_factor} times greater than for the large packages? {is_true}")

solve_candle_arrangement()