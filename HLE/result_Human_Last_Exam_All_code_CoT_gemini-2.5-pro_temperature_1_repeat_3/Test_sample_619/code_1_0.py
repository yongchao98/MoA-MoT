import math

def solve():
    """
    Calculates the minimum value of the sum of set sizes based on the given condition.
    The formula is derived from constructing sets for small n and matching with the given options.
    """
    try:
        n_str = input("Enter the number of sets (n): ")
        n = int(n_str)
        if n < 1:
            print("Please enter a positive integer for n.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The derived formula for the minimum sum is floor(n^2 / 4) + 2.
    n_squared = n * n
    n_squared_by_4 = n_squared / 4
    floor_val = math.floor(n_squared_by_4)
    result = floor_val + 2

    print(f"For n = {n}, the minimum value of the sum is calculated as follows:")
    print(f"Formula: floor(n^2 / 4) + 2")
    print(f"Calculation: floor({n}^2 / 4) + 2 = floor({n_squared} / 4) + 2 = floor({n_squared_by_4}) + 2 = {floor_val} + 2 = {result}")

solve()