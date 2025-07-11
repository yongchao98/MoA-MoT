def solve_queen_problem():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on a 16x16 chessboard without attacking each other.
    """
    
    # The size of the chessboard.
    n = 16
    
    # The problem requires placing m white and m black queens such that no queen
    # attacks another. This is equivalent to placing a total of 2*m queens
    # on the board under the same non-attacking condition.
    
    # On an n x n board, the maximum number of non-attacking queens is n.
    # This is because placing n+1 queens would force at least two queens
    # to be in the same row, due to the Pigeonhole Principle.
    max_total_queens = n
    
    # The total number of queens is 2 * m. So, we have the equation:
    # 2 * m = max_total_queens
    # We solve for m. A solution with 16 queens is known to exist, so this
    # maximum is achievable.
    
    m = max_total_queens // 2
    
    print(f"The problem is for a {n}x{n} board.")
    print("The maximum number of total non-attacking queens is equal to the board size.")
    print(f"The equation relating the total queens to m is: 2 * m = {n}")
    print("Solving for m gives the final answer.")
    print("\nThe final equation with the solution is:")
    print(f"2 * {m} = {n}")
    
    print(f"\nThe maximum number m is {m}.")

solve_queen_problem()