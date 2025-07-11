def solve_queen_problem():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on an N x N chessboard without attacking each other.
    """
    # The size of one side of the chessboard.
    N = 16

    # The maximum number of non-attacking queens on an N x N board is N.
    # This is because if you place more than N queens, by the pigeonhole principle,
    # at least two queens must be in the same row, and would therefore attack each other.
    # It is known that a solution with N queens exists for all N >= 4.
    max_total_queens = N

    # The problem asks for m white and m black queens, for a total of 2*m queens.
    # The condition "without attacking each other" means the entire set of 2*m queens
    # must form a single, non-attacking configuration.
    # Therefore, the total number of queens, 2*m, cannot exceed the maximum possible.
    # 2 * m <= max_total_queens
    # To find the maximum value of m, we use the largest possible number of total queens.
    # 2 * m = max_total_queens
    m = max_total_queens // 2

    # Print the step-by-step reasoning and the final equation.
    print(f"The analysis is for a {N}x{N} chessboard.")
    print(f"The maximum number of total non-attacking queens (K) that can be placed on this board is {max_total_queens}.")
    print("We are placing m white queens and m black queens, making a total of 2*m queens.")
    print("For no queen to attack another, the total number of queens (2*m) must not exceed the maximum K.")
    print("\nTo find the maximum value for m, we form the equation:")
    print(f"2 * m = {max_total_queens}")
    print("\nSolving the equation for m:")
    # The user requested to see the numbers in the final equation.
    print(f"m = {max_total_queens} / 2")
    print(f"m = {m}")
    print(f"\nTherefore, the maximum number m is {m}.")

solve_queen_problem()