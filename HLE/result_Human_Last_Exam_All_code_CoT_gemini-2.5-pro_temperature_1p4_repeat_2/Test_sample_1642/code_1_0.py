def solve_queen_problem():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on a 16x16 chessboard without attacking each other.
    """
    # The size of the chessboard dimension.
    N = 16

    print(f"The analysis is for a {N}x{N} chessboard.")
    print("-" * 30)

    # The problem is to place m white and m black queens without any queen attacking another.
    # This is interpreted as no queen attacking any other queen, regardless of color.
    # Therefore, we have a total of 2*m queens that must all be non-attacking.
    print("The total number of queens to be placed is m white + m black = 2 * m.")

    # The maximum number of non-attacking queens on an N x N board is N (for N>=4).
    max_total_queens = N
    print(f"The maximum number of non-attacking queens on a {N}x{N} board is {max_total_queens}.")

    # To find the maximum value of m, we set the total number of queens (2*m)
    # equal to the maximum possible number of non-attacking queens.
    print("\nThis gives us the final equation to solve for the maximum m:")
    # The prompt requires printing the numbers in the final equation.
    # The equation is 2 * m = 16.
    print("2 * m = 16")

    # Solve the equation for m.
    m = max_total_queens / 2

    print("\nSolving for m:")
    print(f"m = {max_total_queens} / 2")
    print(f"The maximum number m is {int(m)}.")

solve_queen_problem()