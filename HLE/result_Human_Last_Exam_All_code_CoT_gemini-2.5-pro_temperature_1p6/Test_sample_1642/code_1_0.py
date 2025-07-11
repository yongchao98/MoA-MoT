def solve_queen_problem():
    """
    Calculates the maximum number m such that m white queens and m black
    queens can coexist on a 16x16 chessboard without attacking each other.
    """
    # The size of the chessboard.
    N = 16

    # For an N x N board where N >= 4, the maximum number of non-attacking
    # queens, K, that can be placed is equal to N.
    max_total_queens = N

    # The total number of queens on the board is m white queens + m black queens,
    # which is 2 * m. To maximize m, this total must equal the maximum
    # possible number of non-attacking queens.
    # The resulting equation is: 2 * m = max_total_queens
    m = max_total_queens // 2

    print(f"The size of the board is {N}x{N}.")
    print(f"The maximum number of total non-attacking queens on this board is {max_total_queens}.")
    print("This total is composed of m white queens and m black queens (2 * m).")
    print("\nTo find the maximum value of m, we solve the following equation:")
    print(f"2 * m = {max_total_queens}")
    print("\nThe final calculation is:")
    print(f"m = {max_total_queens} / 2 = {m}")
    
    print(f"\nTherefore, the maximum number m is {m}.")

solve_queen_problem()