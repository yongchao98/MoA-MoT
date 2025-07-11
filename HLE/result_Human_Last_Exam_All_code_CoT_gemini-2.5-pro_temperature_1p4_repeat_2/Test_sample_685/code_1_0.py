import sys

def solve_game_probability(n, m):
    """
    Calculates f(n, m), which is 1 if the first player has a winning probability > 50%.
    
    Based on analysis, the winning probability is > 50% if and only if the
    matrix dimensions are not 1x1.
    
    Args:
        n (int): The number of rows in the matrix.
        m (int): The number of columns in the matrix.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        return

    # The condition for the first player having a winning probability > 50%
    # simplifies to checking if the total number of cells is greater than 1.
    total_cells = n * m
    result = 1 if total_cells > 1 else 0
    
    # Print the equation used for the determination
    print(f"To determine f({n}, {m}), we check the condition: {n} * {m} > 1")
    print(f"This evaluates to: {total_cells} > 1, which is {total_cells > 1}.")
    print(f"Therefore, the value of f({n}, {m}) is {result}.")

if __name__ == '__main__':
    # Example usage with n=2, m=3.
    # You can change these values to test other cases.
    n_example = 2
    m_example = 3
    
    # Check if arguments are provided from command line
    if len(sys.argv) == 3:
        try:
            n_example = int(sys.argv[1])
            m_example = int(sys.argv[2])
        except ValueError:
            print("Invalid input. Please provide two integers for n and m.")
            sys.exit(1)

    solve_game_probability(n_example, m_example)