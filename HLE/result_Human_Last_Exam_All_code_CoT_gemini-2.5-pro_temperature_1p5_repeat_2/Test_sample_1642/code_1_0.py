def solve_queen_problem():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on a 16x16 chessboard without attacking each other.
    """
    
    # The size of one side of the chessboard.
    board_size = 16
    
    # The problem states we have 'm' white queens and 'm' black queens.
    # The total number of queens is therefore m + m = 2 * m.
    
    # On an N x N board, you cannot place more than N queens without at least
    # two of them being in the same row, which would cause them to attack each other.
    # Therefore, the maximum total number of non-attacking queens is equal to the board size.
    max_total_queens = board_size
    
    print(f"The chessboard size is {board_size}x{board_size}.")
    print("The total number of queens to be placed is 'm' white + 'm' black = 2*m.")
    print(f"On a {board_size}x{board_size} board, the maximum number of non-attacking queens you can place is {max_total_queens}.")
    print("This is because each queen must occupy a unique row.")
    
    print("\nTo find the maximum value for 'm', we set the total number of queens (2*m) equal to the maximum possible number of queens:")
    
    # We build and print the equation: 2 * m = 16
    # In the equation `2 * m = 16`, the numbers are 2 and 16.
    num_of_colors = 2
    
    print(f"Equation: {num_of_colors} * m = {max_total_queens}")
    
    # Now, we solve for m.
    m = max_total_queens // num_of_colors
    
    print("\nSolving for 'm':")
    print(f"m = {max_total_queens} / {num_of_colors}")
    print(f"m = {m}")
    
    print(f"\nThus, the maximum number m is {m}.")

solve_queen_problem()