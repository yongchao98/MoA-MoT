import sys

def solve():
    """
    This function calculates the maximum value of n for a given k based on the formula n = k^2 - k + 1.
    """
    # Ask the user for the value of k
    try:
        k_str = input("Enter the value of k (integer >= 2): ")
        k = int(k_str)
        if k < 2:
            print("Error: k must be an integer greater than or equal to 2.", file=sys.stderr)
            return
    except ValueError:
        print(f"Error: Invalid input '{k_str}'. Please enter an integer.", file=sys.stderr)
        return

    # Calculate n using the formula n = k^2 - k + 1
    n = k**2 - k + 1
    
    # Print the equation and the result
    print(f"For k = {k}, the maximum value of n is given by the formula:")
    print(f"n = k^2 - k + 1")
    # Output each number in the final equation
    print(f"n = {k}^2 - {k} + 1")
    print(f"n = {k*k} - {k} + 1")
    print(f"n = {n}")

solve()