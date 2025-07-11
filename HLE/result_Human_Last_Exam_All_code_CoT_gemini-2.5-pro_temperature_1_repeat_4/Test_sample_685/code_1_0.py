import sys

def solve_game_probability(n, m):
    """
    This function calculates f(n, m) based on the derived logic.
    f(n, m) = 1 if the first player has a winning position with a probability > 50%.
    f(n, m) = 0 otherwise.

    Args:
        n (int): The number of rows in the matrix.
        m (int): The number of columns in the matrix.
    
    Returns:
        int: The value of f(n, m) which is 0 or 1.
    """
    # The logic boils down to checking if the board is larger than 1x1.
    if n * m > 1:
        # If n*m > 1, P(P-position) < 0.5, so P(win) > 0.5
        result = 1
    else:
        # If n=1 and m=1, P(P-position) = 0.5, so P(win) is not strictly > 0.5
        result = 0
    
    print(f"For a board of size n={n}, m={m}:")
    print(f"The value of f(n, m) is {result}")

def main():
    # Example 1: The base case n=1, m=1
    n1, m1 = 1, 1
    solve_game_probability(n1, m1)
    
    print("-" * 20)
    
    # Example 2: A case where n*m > 1
    n2, m2 = 5, 10
    solve_game_probability(n2, m2)

if __name__ == "__main__":
    main()
