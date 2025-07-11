def solve_queen_problem():
    """
    Calculates the maximum number m such that m white and m black queens
    can coexist on a 16x16 board without attacking each other.
    """
    # Size of the chessboard
    n = 16

    # For an N x N board, the maximum number of non-attacking queens is N.
    # This is because placing N+1 queens would force at least two to be in the same row.
    max_total_queens = n

    # The total number of queens is m white queens + m black queens = 2 * m.
    # We set the total number of queens to the maximum possible value.
    # The equation is: 2 * m = max_total_queens
    # We can now solve for m.
    m = max_total_queens // 2
    
    # We are asked to output the numbers in the final equation.
    # The final equation is 2 * m = 16, which means m = 8.
    # So we print the solved equation.
    
    print(f"On a {n}x{n} chessboard, the maximum number of non-attacking queens is {n}.")
    print("The total number of queens is the sum of white and black queens: m + m = 2 * m.")
    print("To find the maximum value of m, we set up the equation:")
    print(f"2 * m = {max_total_queens}")
    print("\nSolving for m:")
    print(f"{2} * {m} = {max_total_queens}")
    print(f"Therefore, the maximum number m is {m}.")

solve_queen_problem()