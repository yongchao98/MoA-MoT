def solve_queen_problem():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on an N x N chessboard without attacking each other.
    """
    # Define the size of the chessboard
    N = 16

    # The problem involves placing m white queens and m black queens.
    # The total number of queens is 2 * m.
    # The condition "without attacking each other" implies that no queen can attack
    # any other queen on the board, regardless of color.

    # This is equivalent to the classic N-Queens problem, which determines the
    # maximum number of non-attacking queens that can be placed on a board.
    # For an N x N board, the maximum number of non-attacking queens is N.
    max_total_queens = N

    # Therefore, the total number of queens we place (2 * m) must be
    # less than or equal to the maximum possible number of non-attacking queens.
    # This gives us the inequality: 2 * m <= N

    # We solve this inequality for m to find the maximum value.
    # m <= N / 2
    max_m = max_total_queens / 2

    print("Step 1: Understand the constraints.")
    print(f"We have a {N}x{N} chessboard.")
    print("We need to place 'm' white queens and 'm' black queens.")
    print("The condition is that no queen attacks any other queen.")
    print(f"This means the total number of non-attacking queens is m + m = 2 * m.")
    print("-" * 30)

    print("Step 2: Apply the N-Queens problem principle.")
    print(f"The maximum number of non-attacking queens on a {N}x{N} board is {N}.")
    print("-" * 30)

    print("Step 3: Formulate and solve the inequality.")
    print("The total number of queens (2 * m) cannot exceed the maximum possible.")
    print("This gives us the final equation:")
    print(f"2 * m <= {max_total_queens}")
    print(f"m <= {max_total_queens} / 2")
    print(f"m <= {int(max_m)}")
    print("-" * 30)

    print(f"The maximum number m is {int(max_m)}.")

solve_queen_problem()