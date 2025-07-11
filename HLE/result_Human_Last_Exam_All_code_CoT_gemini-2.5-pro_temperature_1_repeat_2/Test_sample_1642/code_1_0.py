def solve_queen_problem():
    """
    Calculates the maximum number m of white/black queens on an nxn board.
    m = k*k, where k is the size of a solid block of queens.
    The condition for a valid placement of two kxk blocks is n - 2k > 2k - 2.
    """
    n = 16
    print(f"The size of the chessboard is {n}x{n}.")
    print("We are looking for the maximum size 'k' of a k x k block of queens.")
    print("The condition for a valid placement is derived from diagonal constraints:")
    print(f"n - 2*k > 2*k - 2")
    print(f"{n} - 2*k > 2*k - 2")
    print("\nLet's solve the inequality for k:")
    print(f"{n} + 2 > 2*k + 2*k")
    print(f"{n + 2} > 4*k")
    print(f"k < {n + 2} / 4")
    print(f"k < {(n + 2) / 4}")

    # The largest integer k must be less than 4.5
    k = (n + 2) // 4
    m = k * k

    print(f"\nThe largest integer value for k is {k}.")
    print(f"The maximum number of queens, m, is k*k.")
    print(f"m = {k} * {k} = {m}")

solve_queen_problem()