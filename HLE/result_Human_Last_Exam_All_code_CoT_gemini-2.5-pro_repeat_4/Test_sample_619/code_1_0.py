import math

def solve():
    """
    Calculates the minimum value of the sum of set sizes based on the derived formula.
    """
    try:
        n_str = input("Enter the value of n (a positive integer): ")
        n = int(n_str)
        if n <= 0:
            print("Please enter a positive integer for n.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The problem asks for the minimum value of sum(|S_i|).
    # Based on constructions for small n, the formula appears to be floor(n^2 / 4) + 2 for n >= 2.
    # For n=1, the minimum sum is 1. We will implement the general formula.
    if n == 1:
        result = 1
        print(f"For n=1, the minimum sum is 1.")
        print(f"The general formula floor(n^2 / 4) + 2 would give floor(1^2 / 4) + 2 = 2.")
    else:
        n_squared = n * n
        val_before_floor = n_squared / 4
        floor_val = math.floor(val_before_floor)
        result = floor_val + 2
        print(f"For n = {n}, the calculation is:")
        print(f"floor({n}^2 / 4) + 2")
        print(f"= floor({n_squared} / 4) + 2")
        print(f"= floor({val_before_floor}) + 2")
        print(f"= {floor_val} + 2")
        print(f"= {result}")

solve()