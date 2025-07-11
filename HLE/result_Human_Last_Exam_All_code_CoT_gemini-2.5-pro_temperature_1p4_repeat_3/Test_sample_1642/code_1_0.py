def solve_queen_problem():
    """
    Calculates the maximum number m of white/black queens on a 16x16 board.
    """
    # Size of the chessboard
    N = 16

    # On an N x N board, the maximum number of non-attacking queens is N.
    max_total_queens = N

    # The total number of queens is m (white) + m (black) = 2 * m.
    # We set this equal to the maximum possible number of queens.
    # Equation: 2 * m = max_total_queens
    # We solve for m.
    m = max_total_queens // 2
    
    # --- Output the reasoning and the final answer ---
    print(f"The size of the chessboard is {N}x{N}.")
    print(f"The maximum total number of non-attacking queens on this board is {max_total_queens}.")
    print("We have 'm' white queens and 'm' black queens, for a total of 2 * m queens.")
    print("To find the maximum value for 'm', we set up the following equation:")
    
    # Print the final equation with each number, as requested.
    num_queen_colors = 2
    print(f"{num_queen_colors} * m = {max_total_queens}")
    
    print("\nSolving for m:")
    print(f"m = {max_total_queens} / {num_queen_colors}")
    print(f"m = {m}")

solve_queen_problem()