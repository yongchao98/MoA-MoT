import math

def solve():
    """
    Calculates the minimum value of the sum of set sizes for a given n.
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

    if n == 1:
        min_sum = 1
        print(f"For n = 1, the condition is trivial for a non-empty set S1.")
        print(f"The minimum size |S1| is 1.")
        print(f"The minimum sum is 1.")
    else:
        # For n >= 2, the formula is floor(n^2 / 4) + 2
        n_squared = n * n
        floor_val = n_squared // 4
        min_sum = floor_val + 2
        
        print(f"For n = {n} (where n >= 2), the formula for the minimum sum is floor(n^2 / 4) + 2.")
        print(f"Calculation:")
        print(f"n^2 = {n}^2 = {n_squared}")
        print(f"n^2 / 4 = {n_squared} / 4 = {n_squared / 4}")
        print(f"floor(n^2 / 4) = floor({n_squared / 4}) = {floor_val}")
        print(f"floor(n^2 / 4) + 2 = {floor_val} + 2 = {min_sum}")
        print(f"\nThe minimum value of the sum is {min_sum}.")

solve()