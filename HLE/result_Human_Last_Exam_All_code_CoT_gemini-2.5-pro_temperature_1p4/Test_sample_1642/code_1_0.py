def solve_queen_problem():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on a 16x16 chessboard without attacking each other.
    """

    # The size of the chessboard.
    board_size = 16

    print(f"The problem is set on a {board_size}x{board_size} chessboard.")
    print("We need to place 'm' white queens and 'm' black queens.")
    print("This makes for a total of 2 * m queens on the board.")
    print("\nThe condition 'without attacking each other' means no queen can attack any other queen, regardless of color.")
    print("This turns the problem into a standard N-Queens puzzle, where N is the board size.")
    
    print(f"\nFor a {board_size}x{board_size} board, the maximum number of non-attacking queens you can place is {board_size}.")
    print("Therefore, the total number of queens (2 * m) must be equal to the board size to be maximum.")

    print("\nWe can set up the following equation:")
    # The variable names 'total_queens' and 'max_queens' are not used in the f-string
    # to ensure the numbers are directly visible in the output per the instructions.
    print(f"2 * m = {board_size}")

    # Solve for m
    m = board_size // 2

    print("\nSolving for m, we get:")
    print(f"m = {board_size} / 2")
    print(f"m = {m}")

    print(f"\nThus, the maximum number m is {m}.")

solve_queen_problem()