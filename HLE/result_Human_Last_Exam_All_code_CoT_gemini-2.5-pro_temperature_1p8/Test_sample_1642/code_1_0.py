def solve_queen_problem():
    """
    Calculates the maximum number m of white and m black queens that can
    coexist on a 16x16 board without attacking each other.
    """
    # The size of the chessboard.
    N = 16

    # For an N x N board where N >= 4, the maximum number of non-attacking queens
    # that can be placed on the board is N. This is a standard result from the
    # n-queens problem.
    max_total_queens = N

    # The problem states we have m white queens and m black queens.
    # The total number of queens is m + m = 2*m.
    # To find the maximum value of m, the total number of queens (2*m) must
    # be equal to the maximum possible non-attacking queens on the board.
    
    # We set up the equation: 2 * m = max_total_queens
    # And solve for m.
    m = max_total_queens / 2

    print(f"The size of the chessboard is N = {N}.")
    print("The condition 'without attacking each other' applies to all queens, regardless of color.")
    print("This means we need to find the maximum number of total non-attacking queens that can be placed on the board.")
    print(f"\nFor a {N}x{N} board, the maximum number of non-attacking queens is {max_total_queens}.")
    print("We are placing m white and m black queens, for a total of 2*m queens.")
    
    print("\nTo find the maximum m, we form the equation:")
    # Printing the equation with its numbers, as requested.
    print(f"2 * m = {max_total_queens}")
    
    print("\nSolving for m:")
    print(f"m = {max_total_queens} / 2")
    
    # The final value of m must be an integer.
    m_solution = int(m)
    print(f"m = {m_solution}")

    print(f"\nThus, the maximum number m is {m_solution}.")

solve_queen_problem()