import sys

def solve():
    """
    Calculates the maximum value of n in terms of k.
    The formula is n = 2k - 1.
    """
    try:
        k_str = input("Enter the value of k (an integer >= 2): ")
        k = int(k_str)
        if k < 2:
            print("Error: k must be an integer greater than or equal to 2.", file=sys.stderr)
            return
    except ValueError:
        print(f"Error: Invalid input '{k_str}'. Please enter an integer.", file=sys.stderr)
        return

    n = 2 * k - 1
    
    print(f"For k = {k}, the maximum value of n is given by the formula:")
    print(f"n_max = 2 * k - 1")
    print(f"n_max = 2 * {k} - 1")
    print(f"n_max = {2*k} - 1")
    print(f"n_max = {n}")

solve()