def solve_queen_problem():
    """
    Calculates the maximum number m for the given queen problem.
    """
    # The size of the chessboard (N x N)
    board_size = 16

    # The number of different colors of queens
    num_colors = 2

    print("Step 1: Understanding the Problem")
    print(f"We are placing 'm' white queens and 'm' black queens on a {board_size}x{board_size} board.")
    print("The condition 'without attacking each other' means no queen can attack any other queen, regardless of color.")
    print("-" * 40)

    print("Step 2: Relating to the N-Queens Problem")
    print("This is equivalent to placing a total of K = m + m = 2*m queens on the board.")
    print("The problem is now to find the maximum number of non-attacking queens on a 16x16 board.")
    print("-" * 40)

    print("Step 3: Applying the Fundamental Limit")
    print(f"For an {board_size}x{board_size} board, the maximum number of non-attacking queens you can place is {board_size}.")
    print("This is because each queen must be in a unique row, and there are only 16 rows.")
    max_total_queens = board_size
    print(f"Therefore, the total number of queens (2 * m) cannot exceed {max_total_queens}.")
    print("-" * 40)

    print("Step 4: Solving for m")
    # The equation is: num_colors * m = max_total_queens
    m = max_total_queens // num_colors
    print(f"We have the equation: {num_colors} * m = {max_total_queens}")
    print(f"Solving for m gives: m = {max_total_queens} / {num_colors}")
    print(f"The maximum value for m is {m}.")
    print("-" * 40)

    print("Final Equation:")
    # The final output must show the equation with the numbers plugged in.
    print(f"{num_colors} * {m} = {max_total_queens}")

solve_queen_problem()