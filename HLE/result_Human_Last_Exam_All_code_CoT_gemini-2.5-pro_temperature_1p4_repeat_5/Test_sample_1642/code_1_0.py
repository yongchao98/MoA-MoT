def solve_queen_problem():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on a 16x16 chessboard without attacking each other.
    """
    # Step 1: Define the size of the chessboard.
    n = 16

    # Step 2: The maximum number of non-attacking queens on an n x n board is n.
    total_queens = n

    # Step 3: The total number of queens is m white + m black, so total = 2 * m.
    # We need to solve the equation: 2 * m = n
    m = total_queens // 2

    # Step 4: Print the reasoning and the final equation.
    print(f"The problem is to place 'm' white and 'm' black queens on a {n}x{n} board with no attacks.")
    print("This is equivalent to the N-Queens problem, where the maximum number of non-attacking queens on an NxN board is N.")
    print(f"For a {n}x{n} board, the maximum total number of non-attacking queens is {total_queens}.")
    print("\nSince the total number of queens is m + m, we have the following equation:")
    print(f"m + m = {total_queens}")
    print(f"2 * m = {total_queens}")
    print(f"m = {total_queens} / 2")
    print(f"m = {m}")
    print(f"\nThus, the maximum number m is {m}.")

solve_queen_problem()