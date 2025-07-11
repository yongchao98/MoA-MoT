import math

def solve():
    """
    Calculates the minimum value of the sum of set sizes based on the derived formula.
    """
    try:
        n_str = input("Enter the number of sets (n): ")
        n = int(n_str)
        if n <= 0:
            print("Please enter a positive integer for n.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The minimum value is given by the formula floor(n^2 / 4) + 2.
    val_floor = math.floor(n**2 / 4)
    min_sum = val_floor + 2

    print(f"For n = {n}, the minimum value of the sum is:")
    # The prompt asks to "output each number in the final equation"
    print(f"floor({n}^2 / 4) + 2 = {val_floor} + 2 = {min_sum}")

solve()