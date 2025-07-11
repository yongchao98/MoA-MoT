def solve_queen_problem():
    """
    Calculates the maximum number m such that m white and m black queens
    can coexist on a 16x16 board without attacking each other.
    """
    # The size of the chessboard side.
    N = 16

    # The problem is to place m white and m black queens such that no queen
    # attacks another. This is a total of 2*m queens.
    # The maximum number of non-attacking queens on an N x N board is N.
    max_total_queens = N

    # This gives us the inequality: 2 * m <= N
    # To find the maximum m, we solve for m.
    max_m = max_total_queens // 2

    print(f"The board size is {N}x{N}.")
    print("The goal is to place m white and m black queens without any queen attacking another.")
    print(f"This means the total number of non-attacking queens is 2 * m.")
    print("\n")
    print(f"The maximum number of non-attacking queens on a {N}x{N} board is {max_total_queens}.")
    print(f"Therefore, the total number of queens must satisfy: 2 * m <= {max_total_queens}.")
    print("\n")
    print("To find the maximum value for m, we solve the equation derived from the limit:")
    # As requested, printing each number in the final equation
    print(f"m = {max_total_queens} / 2")
    print(f"m = {max_m}")
    print("\n")
    print(f"The maximum number m is {max_m}.")

solve_queen_problem()